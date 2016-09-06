/*
 * my_pv_system_3_ss_controller.h
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
#ifndef RTW_HEADER_my_pv_system_3_ss_controller_h_
#define RTW_HEADER_my_pv_system_3_ss_controller_h_
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <math.h>
#ifndef my_pv_system_3_ss_controller_COMMON_INCLUDES_
# define my_pv_system_3_ss_controller_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* my_pv_system_3_ss_controller_COMMON_INCLUDES_ */

#include "my_pv_system_3_ss_controller_types.h"

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
#define my_pv_system_3_ss_controller_rtModel RT_MODEL_my_pv_system_3_ss_controller_T

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
  creal_T RealImagtoComplex;           /* '<S16>/Real-Imag to Complex' */
  creal_T MagnitudeAngletoComplex;     /* '<S11>/Magnitude-Angle to Complex' */
  creal_T RealImagtoComplex_n;         /* '<S13>/Real-Imag to Complex' */
  creal_T RealImagtoComplex_l;         /* '<S36>/Real-Imag to Complex' */
  real_T SFunction;                    /* '<S1>/S-Function' */
  real_T Sum;                          /* '<S1>/Sum' */
  real_T Gate[4];                      /* '<S2>/Memory' */
  real_T Iph;                          /* '<S2>/Memory' */
  real_T Io;                           /* '<S2>/Memory' */
  real_T OpMonitor_o1;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o2;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o3;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o4;                 /* '<S2>/OpMonitor' */
  real_T SFunction_h[11];              /* '<S91>/S-Function' */
  real_T RateTransition1;              /* '<S2>/Rate Transition1' */
  real_T Apu;                          /* '<S11>/A->pu' */
  real_T UnitDelay;                    /* '<S31>/Unit Delay' */
  real_T avoiddivisionbyzero;          /* '<S20>/avoid division by zero' */
  real_T MathFunction;                 /* '<S20>/Math Function' */
  real_T Gain;                         /* '<S20>/Gain' */
  real_T SFunction_o;                  /* '<S57>/S-Function' */
  real_T DiscreteTimeIntegrator;       /* '<S31>/Discrete-Time Integrator' */
  real_T MathFunction_i;               /* '<S31>/Math Function' */
  real_T FirstcycleofsimulationId092Iq0;/* '<S20>/First cycle of simulation Id=0.92, Iq=0' */
  real_T Switch[2];                    /* '<S20>/Switch' */
  real_T Duk[2];                       /* '<S52>/D*u(k)' */
  real_T x1k[2];                       /* '<S52>/Delay_x1' */
  real_T C11[2];                       /* '<S55>/C11' */
  real_T x2k[2];                       /* '<S52>/Delay_x2' */
  real_T C12[2];                       /* '<S55>/C12' */
  real_T sum2[2];                      /* '<S55>/sum2' */
  real_T yk[2];                        /* '<S52>/C*X(k)+D*u(k)' */
  real_T UnitDelay2;                   /* '<S4>/Unit Delay2' */
  real_T Sum_e[2];                     /* '<S8>/Sum' */
  real_T ProportionalGain[2];          /* '<S15>/Proportional Gain' */
  real_T Integrator[2];                /* '<S15>/Integrator' */
  real_T Sum_g[2];                     /* '<S15>/Sum' */
  real_T Saturate[2];                  /* '<S15>/Saturate' */
  real_T RateTransition;               /* '<S2>/Rate Transition' */
  real_T Vpu;                          /* '<S11>/V->pu' */
  real_T TrigonometricFunction;        /* '<S16>/Trigonometric Function' */
  real_T Gain1;                        /* '<S16>/Gain1' */
  real_T Product1;                     /* '<S16>/Product1' */
  real_T Integ4;                       /* '<S23>/Integ4' */
  real_T Freq;                         /* '<S23>/To avoid division  by zero' */
  real_T Numberofsamplespercycle;      /* '<S23>/Number of samples per cycle' */
  real_T RoundingFunction;             /* '<S23>/Rounding Function' */
  real_T Delay;                        /* '<S23>/Gain' */
  real_T SFunction_p;                  /* '<S25>/S-Function' */
  real_T UnitDelay_f;                  /* '<S24>/Unit Delay' */
  real_T DigitalClock;                 /* '<S23>/Digital  Clock' */
  real_T UnitDelay1;                   /* '<S23>/Unit Delay1' */
  real_T Switch_g;                     /* '<S23>/Switch' */
  real_T TrigonometricFunction3;       /* '<S16>/Trigonometric Function3' */
  real_T Gain3;                        /* '<S16>/Gain3' */
  real_T Product2;                     /* '<S16>/Product2' */
  real_T Integ4_h;                     /* '<S26>/Integ4' */
  real_T Freq_m;                       /* '<S26>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_m;    /* '<S26>/Number of samples per cycle' */
  real_T RoundingFunction_p;           /* '<S26>/Rounding Function' */
  real_T Delay_d;                      /* '<S26>/Gain' */
  real_T SFunction_i;                  /* '<S28>/S-Function' */
  real_T UnitDelay_l;                  /* '<S27>/Unit Delay' */
  real_T DigitalClock_m;               /* '<S26>/Digital  Clock' */
  real_T UnitDelay1_e;                 /* '<S26>/Unit Delay1' */
  real_T Switch_k;                     /* '<S26>/Switch' */
  real_T ComplextoMagnitudeAngle_o1;   /* '<S16>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2;   /* '<S16>/Complex to Magnitude-Angle' */
  real_T RadDeg;                       /* '<S16>/Rad->Deg.' */
  real_T torad;                        /* '<S11>/to-rad' */
  real_T ComplextoRealImag_o1;         /* '<S11>/Complex to Real-Imag' */
  real_T ComplextoRealImag_o2;         /* '<S11>/Complex to Real-Imag' */
  real_T Rff;                          /* '<S8>/Rff ' */
  real_T Lff;                          /* '<S8>/Lff  ' */
  real_T Feedforward;                  /* '<S8>/Add1' */
  real_T Rff_e;                        /* '<S8>/Rff' */
  real_T Lff_i;                        /* '<S8>/Lff' */
  real_T Add3;                         /* '<S8>/Add3' */
  real_T Add2[2];                      /* '<S8>/Add2' */
  real_T IntegralGain[2];              /* '<S15>/Integral Gain' */
  real_T Saturation[2];                /* '<S8>/Saturation' */
  real_T RateTransition3;              /* '<S2>/Rate Transition3' */
  real_T RateTransition4;              /* '<S2>/Rate Transition4' */
  real_T DigitalClock_a;               /* '<S29>/Digital  Clock' */
  real_T RateTransition2;              /* '<S2>/Rate Transition2' */
  real_T Integ4_c;                     /* '<S29>/Integ4' */
  real_T K1;                           /* '<S29>/K1' */
  real_T SFunction_pd;                 /* '<S30>/S-Function' */
  real_T UnitDelay_h;                  /* '<S29>/Unit Delay' */
  real_T UnitDelay1_p;                 /* '<S29>/Unit Delay1' */
  real_T Switch_c;                     /* '<S29>/Switch' */
  real_T TrigonometricFunction2;       /* '<S31>/Trigonometric Function2' */
  real_T Product1_c;                   /* '<S31>/Product1' */
  real_T Integ4_e;                     /* '<S45>/Integ4' */
  real_T Freq_c;                       /* '<S45>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_f;    /* '<S45>/Number of samples per cycle' */
  real_T RoundingFunction_b;           /* '<S45>/Rounding Function' */
  real_T Delay_m;                      /* '<S45>/Gain' */
  real_T SFunction_e;                  /* '<S47>/S-Function' */
  real_T UnitDelay_o;                  /* '<S46>/Unit Delay' */
  real_T DigitalClock_e;               /* '<S45>/Digital  Clock' */
  real_T UnitDelay1_a;                 /* '<S45>/Unit Delay1' */
  real_T Switch_d;                     /* '<S45>/Switch' */
  real_T Divide;                       /* '<S31>/Divide' */
  real_T DiscreteDerivative;           /* '<S33>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_g;     /* '<S33>/Discrete-Time Integrator' */
  real_T Kp4;                          /* '<S33>/Kp4' */
  real_T Sum6;                         /* '<S33>/Sum6' */
  real_T Saturation1;                  /* '<S33>/Saturation1' */
  real_T Gain10;                       /* '<S31>/Gain10' */
  real_T RateLimiter;                  /* '<S31>/Rate Limiter' */
  real_T x1k_d;                        /* '<S48>/Delay_x1' */
  real_T A11;                          /* '<S49>/A11' */
  real_T x2k_k;                        /* '<S48>/Delay_x2' */
  real_T A12;                          /* '<S49>/A12' */
  real_T A21;                          /* '<S49>/A21' */
  real_T A22;                          /* '<S49>/A22' */
  real_T sum2_c;                       /* '<S49>/sum2' */
  real_T sum3;                         /* '<S49>/sum3' */
  real_T B11;                          /* '<S50>/B11' */
  real_T x1k1;                         /* '<S48>/A*x1(k) + B*u1(k) ' */
  real_T B21;                          /* '<S50>/B21' */
  real_T x2k1;                         /* '<S48>/A*x2(k) + B*u2(k)' */
  real_T Duk_b;                        /* '<S48>/D*u(k)' */
  real_T C11_d;                        /* '<S51>/C11' */
  real_T C12_g;                        /* '<S51>/C12' */
  real_T sum2_o;                       /* '<S51>/sum2' */
  real_T yk_e;                         /* '<S48>/C*X(k)+D*u(k)' */
  real_T A11_e[2];                     /* '<S53>/A11' */
  real_T A12_i[2];                     /* '<S53>/A12' */
  real_T A21_o[2];                     /* '<S53>/A21' */
  real_T A22_p[2];                     /* '<S53>/A22' */
  real_T sum2_e[2];                    /* '<S53>/sum2' */
  real_T sum3_g[2];                    /* '<S53>/sum3' */
  real_T B11_o[2];                     /* '<S54>/B11' */
  real_T x1k1_l[2];                    /* '<S52>/A*x1(k) + B*u1(k) ' */
  real_T B21_p[2];                     /* '<S54>/B21' */
  real_T x2k1_k[2];                    /* '<S52>/A*x2(k) + B*u2(k)' */
  real_T Constant1;                    /* '<S20>/Constant1' */
  real_T Add3_n;                       /* '<S62>/Add3' */
  real_T DigitalClock_mz;              /* '<S84>/Digital Clock' */
  real_T Add1;                         /* '<S84>/Add1' */
  real_T MathFunction_b;               /* '<S84>/Math Function' */
  real_T ib1;                          /* '<S84>/1\ib1' */
  real_T LookupTable;                  /* '<S84>/Lookup Table' */
  real_T Add3_o;                       /* '<S84>/Add3' */
  real_T Gain1_p;                      /* '<S62>/Gain1' */
  real_T MUL1;                         /* '<S62>/MUL1' */
  real_T Add4;                         /* '<S62>/Add4' */
  real_T UnitDelay_m;                  /* '<S4>/Unit Delay' */
  real_T Gain_a;                       /* '<S67>/Gain' */
  real_T DataTypeConversion[4];        /* '<S12>/Data Type Conversion' */
  real_T Switch_d2;                    /* '<S4>/Switch' */
  real_T Add1_a;                       /* '<S13>/Add1' */
  real_T UnitDelay3[2];                /* '<S4>/Unit Delay3' */
  real_T Gain1_g;                      /* '<S13>/Gain1' */
  real_T Product;                      /* '<S13>/Product' */
  real_T Product1_d[2];                /* '<S13>/Product1' */
  real_T ComplextoMagnitudeAngle_o1_n; /* '<S13>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_p; /* '<S13>/Complex to Magnitude-Angle' */
  real_T Add2_f;                       /* '<S13>/Add2' */
  real_T TrigonometricFunction_h;      /* '<S13>/Trigonometric Function' */
  real_T Product2_n;                   /* '<S13>/Product2' */
  real_T UnitDelay1_c;                 /* '<S4>/Unit Delay1' */
  real_T Sum_k;                        /* '<S14>/Sum' */
  real_T Rtot_pu2;                     /* '<S14>/Rtot_pu2' */
  real_T IntegralGain_a;               /* '<S86>/Integral Gain' */
  real_T Integrator_e;                 /* '<S86>/Integrator' */
  real_T ProportionalGain_b;           /* '<S86>/Proportional Gain' */
  real_T Sum_i;                        /* '<S86>/Sum' */
  real_T Saturate_p;                   /* '<S86>/Saturate' */
  real_T Iph_p;                        /* '<S2>/MATLAB Function' */
  real_T Io_c;                         /* '<S2>/MATLAB Function' */
  real_T Fcn;                          /* '<S61>/Fcn' */
  real_T Fcn1;                         /* '<S61>/Fcn1' */
  real_T Fcn_l;                        /* '<S60>/Fcn' */
  real_T Fcn1_f;                       /* '<S60>/Fcn1' */
  real_T Switch_p[2];                  /* '<S56>/Switch' */
  real_T Sum1;                         /* '<S46>/Sum1' */
  real_T Sum5;                         /* '<S46>/Sum5' */
  real_T Product5;                     /* '<S46>/Product5' */
  real_T Gain1_pw;                     /* '<S46>/Gain1' */
  real_T Sum4;                         /* '<S46>/Sum4' */
  real_T Product2_b;                   /* '<S46>/Product2' */
  real_T Product4;                     /* '<S46>/Product4' */
  real_T Sum7;                         /* '<S45>/Sum7' */
  real_T Meanvalue;                    /* '<S45>/Product' */
  real_T Sum5_b;                       /* '<S45>/Sum5' */
  real_T TrigonometricFunction_p;      /* '<S36>/Trigonometric Function' */
  real_T Gain1_i;                      /* '<S36>/Gain1' */
  real_T Product1_o;                   /* '<S36>/Product1' */
  real_T Integ4_d;                     /* '<S39>/Integ4' */
  real_T Freq_h;                       /* '<S39>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_h;    /* '<S39>/Number of samples per cycle' */
  real_T RoundingFunction_e;           /* '<S39>/Rounding Function' */
  real_T Delay_db;                     /* '<S39>/Gain' */
  real_T SFunction_ij;                 /* '<S41>/S-Function' */
  real_T UnitDelay_b;                  /* '<S40>/Unit Delay' */
  real_T DigitalClock_o;               /* '<S39>/Digital  Clock' */
  real_T UnitDelay1_c2;                /* '<S39>/Unit Delay1' */
  real_T Switch_j;                     /* '<S39>/Switch' */
  real_T TrigonometricFunction3_p;     /* '<S36>/Trigonometric Function3' */
  real_T Gain3_c;                      /* '<S36>/Gain3' */
  real_T Product2_i;                   /* '<S36>/Product2' */
  real_T Integ4_dy;                    /* '<S42>/Integ4' */
  real_T Freq_mz;                      /* '<S42>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_d;    /* '<S42>/Number of samples per cycle' */
  real_T RoundingFunction_h;           /* '<S42>/Rounding Function' */
  real_T Delay_i;                      /* '<S42>/Gain' */
  real_T SFunction_ed;                 /* '<S44>/S-Function' */
  real_T UnitDelay_fb;                 /* '<S43>/Unit Delay' */
  real_T DigitalClock_l;               /* '<S42>/Digital  Clock' */
  real_T UnitDelay1_j;                 /* '<S42>/Unit Delay1' */
  real_T Switch_m;                     /* '<S42>/Switch' */
  real_T ComplextoMagnitudeAngle_o1_d; /* '<S36>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_i; /* '<S36>/Complex to Magnitude-Angle' */
  real_T RadDeg_h;                     /* '<S36>/Rad->Deg.' */
  real_T Saturation_n;                 /* '<S32>/Saturation' */
  real_T MathFunction_k;               /* '<S32>/Math Function' */
  real_T Sum1_c;                       /* '<S43>/Sum1' */
  real_T Sum5_o;                       /* '<S43>/Sum5' */
  real_T Product5_e;                   /* '<S43>/Product5' */
  real_T Gain1_n;                      /* '<S43>/Gain1' */
  real_T Sum4_g;                       /* '<S43>/Sum4' */
  real_T Product2_m;                   /* '<S43>/Product2' */
  real_T Product4_d;                   /* '<S43>/Product4' */
  real_T Sum7_m;                       /* '<S42>/Sum7' */
  real_T Meanvalue_g;                  /* '<S42>/Product' */
  real_T Sum5_i;                       /* '<S42>/Sum5' */
  real_T Sum1_p;                       /* '<S40>/Sum1' */
  real_T Sum5_k;                       /* '<S40>/Sum5' */
  real_T Product5_n;                   /* '<S40>/Product5' */
  real_T Gain1_o;                      /* '<S40>/Gain1' */
  real_T Sum4_n;                       /* '<S40>/Sum4' */
  real_T Product2_d;                   /* '<S40>/Product2' */
  real_T Product4_o;                   /* '<S40>/Product4' */
  real_T Sum7_e;                       /* '<S39>/Sum7' */
  real_T Meanvalue_h;                  /* '<S39>/Product' */
  real_T Sum5_f;                       /* '<S39>/Sum5' */
  real_T Gain1_c;                      /* '<S29>/Gain1' */
  real_T Gain_f;                       /* '<S29>/Gain' */
  real_T Correction;                   /* '<S29>/Sum1' */
  real_T Sum7_i;                       /* '<S29>/Sum7' */
  real_T Mean;                         /* '<S29>/Product' */
  real_T Sum5_c;                       /* '<S29>/Sum5' */
  real_T Sum1_cr;                      /* '<S27>/Sum1' */
  real_T Sum5_g;                       /* '<S27>/Sum5' */
  real_T Product5_a;                   /* '<S27>/Product5' */
  real_T Gain1_if;                     /* '<S27>/Gain1' */
  real_T Sum4_b;                       /* '<S27>/Sum4' */
  real_T Product2_dp;                  /* '<S27>/Product2' */
  real_T Product4_g;                   /* '<S27>/Product4' */
  real_T Sum7_i1;                      /* '<S26>/Sum7' */
  real_T Meanvalue_e;                  /* '<S26>/Product' */
  real_T Sum5_g2;                      /* '<S26>/Sum5' */
  real_T Sum1_c3;                      /* '<S24>/Sum1' */
  real_T Sum5_cr;                      /* '<S24>/Sum5' */
  real_T Product5_j;                   /* '<S24>/Product5' */
  real_T Gain1_b;                      /* '<S24>/Gain1' */
  real_T Sum4_e;                       /* '<S24>/Sum4' */
  real_T Product2_f;                   /* '<S24>/Product2' */
  real_T Product4_c;                   /* '<S24>/Product4' */
  real_T Sum7_eg;                      /* '<S23>/Sum7' */
  real_T Meanvalue_c;                  /* '<S23>/Product' */
  real_T Sum5_m;                       /* '<S23>/Sum5' */
  real_T TmpSignalConversionAtSFunctionI[4];/* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T D;                            /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  uint8_T Compare;                     /* '<S58>/Compare' */
  uint8_T Compare_e;                   /* '<S59>/Compare' */
  boolean_T RelationalOperator;        /* '<S23>/Relational Operator' */
  boolean_T RelationalOperator_m;      /* '<S26>/Relational Operator' */
  boolean_T RelationalOperator_j;      /* '<S29>/Relational Operator' */
  boolean_T RelationalOperator_l;      /* '<S45>/Relational Operator' */
  boolean_T RelationalOperator1;       /* '<S67>/Relational Operator1' */
  boolean_T RelationalOperator3;       /* '<S67>/Relational Operator3' */
  boolean_T LogicalOperator4[2];       /* '<S12>/Logical Operator4' */
  boolean_T RelationalOperator_d;      /* '<S39>/Relational Operator' */
  boolean_T RelationalOperator_k;      /* '<S42>/Relational Operator' */
} B_my_pv_system_3_ss_controller_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<S31>/Unit Delay' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S31>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE[2];           /* '<S52>/Delay_x1' */
  real_T Delay_x2_DSTATE[2];           /* '<S52>/Delay_x2' */
  real_T UnitDelay2_DSTATE;            /* '<S4>/Unit Delay2' */
  real_T Integrator_DSTATE[2];         /* '<S15>/Integrator' */
  real_T Integ4_DSTATE;                /* '<S23>/Integ4' */
  real_T UnitDelay_DSTATE_g;           /* '<S24>/Unit Delay' */
  real_T UnitDelay1_DSTATE;            /* '<S23>/Unit Delay1' */
  real_T Integ4_DSTATE_d;              /* '<S26>/Integ4' */
  real_T UnitDelay_DSTATE_i;           /* '<S27>/Unit Delay' */
  real_T UnitDelay1_DSTATE_b;          /* '<S26>/Unit Delay1' */
  real_T Integ4_DSTATE_d5;             /* '<S29>/Integ4' */
  real_T UnitDelay_DSTATE_iq;          /* '<S29>/Unit Delay' */
  real_T UnitDelay1_DSTATE_a;          /* '<S29>/Unit Delay1' */
  real_T Integ4_DSTATE_k;              /* '<S45>/Integ4' */
  real_T UnitDelay_DSTATE_a;           /* '<S46>/Unit Delay' */
  real_T UnitDelay1_DSTATE_i;          /* '<S45>/Unit Delay1' */
  real_T DiscreteDerivative_states;    /* '<S33>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_DSTATE_d;/* '<S33>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE_d;            /* '<S48>/Delay_x1' */
  real_T Delay_x2_DSTATE_m;            /* '<S48>/Delay_x2' */
  real_T UnitDelay_DSTATE_il;          /* '<S4>/Unit Delay' */
  real_T UnitDelay3_DSTATE[2];         /* '<S4>/Unit Delay3' */
  real_T UnitDelay1_DSTATE_e;          /* '<S4>/Unit Delay1' */
  real_T Integrator_DSTATE_d;          /* '<S86>/Integrator' */
  real_T Integ4_DSTATE_ky;             /* '<S39>/Integ4' */
  real_T UnitDelay_DSTATE_l;           /* '<S40>/Unit Delay' */
  real_T UnitDelay1_DSTATE_o;          /* '<S39>/Unit Delay1' */
  real_T Integ4_DSTATE_j;              /* '<S42>/Integ4' */
  real_T UnitDelay_DSTATE_d;           /* '<S43>/Unit Delay' */
  real_T UnitDelay1_DSTATE_k;          /* '<S42>/Unit Delay1' */
  real_T SFunction_PreviousInput;      /* '<S1>/S-Function' */
  real_T Memory_1_PreviousInput[4];    /* '<S2>/Memory' */
  real_T Memory_2_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_3_PreviousInput;       /* '<S2>/Memory' */
  real_T DiscreteDerivative_tmp;       /* '<S33>/Discrete Derivative ' */
  real_T PrevY;                        /* '<S31>/Rate Limiter' */
  real_T Vold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Pold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Dold;                         /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T SFunction_RWORK;              /* '<S57>/S-Function' */
  real_T SFunction_RWORK_a;            /* '<S25>/S-Function' */
  real_T SFunction_RWORK_o;            /* '<S28>/S-Function' */
  real_T SFunction_RWORK_k;            /* '<S30>/S-Function' */
  real_T SFunction_RWORK_j;            /* '<S47>/S-Function' */
  real_T SFunction_RWORK_kx;           /* '<S41>/S-Function' */
  real_T SFunction_RWORK_m;            /* '<S44>/S-Function' */
  void *OpMonitor_PWORK;               /* '<S2>/OpMonitor' */
  void *SFunction_PWORK;               /* '<S91>/S-Function' */
  void *SFunction_PWORK_a;             /* '<S57>/S-Function' */
  struct {
    void *LoggedData;
  } PI_Ireg1_PWORK;                    /* '<S8>/PI_Ireg1' */

  void *SFunction_PWORK_c;             /* '<S25>/S-Function' */
  void *SFunction_PWORK_f;             /* '<S28>/S-Function' */
  void *SFunction_PWORK_i;             /* '<S30>/S-Function' */
  void *SFunction_PWORK_a5;            /* '<S47>/S-Function' */
  void *SFunction_PWORK_n;             /* '<S41>/S-Function' */
  void *SFunction_PWORK_p;             /* '<S44>/S-Function' */
  int_T SFunction_IWORK[5];            /* '<S89>/S-Function' */
  int_T SFunction_IWORK_p;             /* '<S57>/S-Function' */
  int_T SFunction_IWORK_f;             /* '<S25>/S-Function' */
  int_T SFunction_IWORK_pq;            /* '<S28>/S-Function' */
  int_T SFunction_IWORK_b;             /* '<S30>/S-Function' */
  int_T SFunction_IWORK_c;             /* '<S47>/S-Function' */
  int_T SFunction_IWORK_i[5];          /* '<S90>/S-Function' */
  int_T SFunction_IWORK_pd;            /* '<S41>/S-Function' */
  int_T SFunction_IWORK_c4;            /* '<S44>/S-Function' */
  uint8_T Integ4_SYSTEM_ENABLE;        /* '<S23>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_e;      /* '<S26>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_d;      /* '<S29>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_k;      /* '<S45>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_g;      /* '<S39>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_gx;     /* '<S42>/Integ4' */
  boolean_T Vold_not_empty;            /* '<S4>/MPPT Controller using Perturbe  & Observe technique  ' */
  boolean_T AutomaticGainControl_MODE; /* '<S31>/Automatic Gain Control' */
} DW_my_pv_system_3_ss_controller_T;

/* Backward compatible GRT Identifiers */
#define rtB                            my_pv_system_3_ss_controller_B
#define BlockIO                        B_my_pv_system_3_ss_controller_T
#define rtP                            my_pv_system_3_ss_controller_P
#define Parameters                     P_my_pv_system_3_ss_controller_T
#define rtDWork                        my_pv_system_3_ss_controller_DW
#define D_Work                         DW_my_pv_system_3_ss_controller_T

/* Parameters (auto storage) */
struct P_my_pv_system_3_ss_controller_T_ {
  real_T PLL_AGC;                      /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S31>/Constant1'
                                        */
  real_T AlphaBetaZerotodq0_Alignment; /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S56>/Constant'
                                        */
  real_T InverterControl_Fnom;         /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S13>/Constant4'
                                        *   '<S20>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  real_T InverterControl_Increment_MPPT;/* Mask Parameter: InverterControl_Increment_MPPT
                                         * Referenced by: '<S9>/Iph_3'
                                         */
  real_T Discrete_Init;                /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S33>/Discrete-Time Integrator'
                                        */
  real_T Discrete_Kd;                  /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S33>/Discrete Derivative '
                                        */
  real_T InverterControl_Ki_Ireg;      /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S15>/Integral Gain'
                                        */
  real_T InverterControl_Ki_VDCreg;    /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S86>/Integral Gain'
                                        */
  real_T Discrete_Kp;                  /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S33>/Kp4'
                                        */
  real_T InverterControl_Kp_Ireg;      /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S15>/Proportional Gain'
                                        */
  real_T InverterControl_Kp_VDCreg;    /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S86>/Proportional Gain'
                                        */
  real_T PI_LowerSaturationLimit;      /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S15>/Saturate'
                                        */
  real_T PI_LowerSaturationLimit_p;    /* Mask Parameter: PI_LowerSaturationLimit_p
                                        * Referenced by: '<S86>/Saturate'
                                        */
  real_T PWM_Generator_MinMax[2];      /* Mask Parameter: PWM_Generator_MinMax
                                        * Referenced by: '<S12>/Constant10'
                                        */
  real_T InverterControl_Pnom;         /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S11>/A->pu'
                                        */
  real_T InverterControl_Ts_Control;   /* Mask Parameter: InverterControl_Ts_Control
                                        * Referenced by:
                                        *   '<S13>/Constant4'
                                        *   '<S23>/Gain'
                                        *   '<S26>/Gain'
                                        *   '<S45>/Gain'
                                        *   '<S39>/Gain'
                                        *   '<S42>/Gain'
                                        */
  real_T PI_UpperSaturationLimit;      /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S15>/Saturate'
                                        */
  real_T PI_UpperSaturationLimit_f;    /* Mask Parameter: PI_UpperSaturationLimit_f
                                        * Referenced by: '<S86>/Saturate'
                                        */
  real_T InverterControl_Vdc_ref_Init; /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S4>/Vnom_dc1'
                                        *   '<S9>/Iph_'
                                        */
  real_T InverterControl_Vnom_dc;      /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S14>/Rtot_pu2'
                                        */
  real_T InverterControl_Vnom_prim;    /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S11>/A->pu'
                                        *   '<S11>/V->pu'
                                        *   '<S13>/Constant3'
                                        */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S58>/Constant'
                                        */
  real_T CompareToConstant1_const;     /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S59>/Constant'
                                        */
  real_T Gain1_Gain;                   /* Expression: 0.5
                                        * Referenced by: '<S24>/Gain1'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 0.5
                                        * Referenced by: '<S27>/Gain1'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: sps.K2
                                        * Referenced by: '<S29>/Gain1'
                                        */
  real_T Gain_Gain;                    /* Expression: sps.K1
                                        * Referenced by: '<S29>/Gain'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: 0.5
                                        * Referenced by: '<S40>/Gain1'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 0.5
                                        * Referenced by: '<S43>/Gain1'
                                        */
  real_T Gain_Y0;                      /* Expression: [1]
                                        * Referenced by: '<S32>/Gain'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 2
                                        * Referenced by: '<S36>/Gain1'
                                        */
  real_T Integ4_gainval;               /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S39>/Integ4'
                                        */
  real_T Integ4_IC;                    /* Expression: 0
                                        * Referenced by: '<S39>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSat;/* Expression: 1e6
                                         * Referenced by: '<S39>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_LowerSat;/* Expression: eps
                                         * Referenced by: '<S39>/To avoid division  by zero'
                                         */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: MaxDelay
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: Ts
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: InitialValue
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P4_Size[2];         /* Computed Parameter: SFunction_P4_Size
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T SFunction_P4;                 /* Expression: DFT
                                        * Referenced by: '<S41>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S40>/Unit Delay'
                                        */
  real_T Constant_Value;               /* Expression: 1/sps.Finit
                                        * Referenced by: '<S39>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: sps.Vinit
                                        * Referenced by: '<S39>/Unit Delay1'
                                        */
  real_T Gain3_Gain;                   /* Expression: 2
                                        * Referenced by: '<S36>/Gain3'
                                        */
  real_T Integ4_gainval_h;             /* Computed Parameter: Integ4_gainval_h
                                        * Referenced by: '<S42>/Integ4'
                                        */
  real_T Integ4_IC_k;                  /* Expression: 0
                                        * Referenced by: '<S42>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_n;/* Expression: 1e6
                                          * Referenced by: '<S42>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_m;/* Expression: eps
                                          * Referenced by: '<S42>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_a[2];       /* Computed Parameter: SFunction_P1_Size_a
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P1_a;               /* Expression: MaxDelay
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P2_Size_k[2];       /* Computed Parameter: SFunction_P2_Size_k
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P2_l;               /* Expression: Ts
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P3_Size_h[2];       /* Computed Parameter: SFunction_P3_Size_h
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P3_f;               /* Expression: InitialValue
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P4_Size_l[2];       /* Computed Parameter: SFunction_P4_Size_l
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P4_k;               /* Expression: DFT
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_c; /* Expression: 0
                                        * Referenced by: '<S43>/Unit Delay'
                                        */
  real_T Constant_Value_c;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S42>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_c;/* Expression: sps.Vinit
                                        * Referenced by: '<S42>/Unit Delay1'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S36>/Rad->Deg.'
                                        */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S32>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: eps
                                        * Referenced by: '<S32>/Saturation'
                                        */
  real_T Gain1_Gain_hn;                /* Expression: 0.5
                                        * Referenced by: '<S46>/Gain1'
                                        */
  real_T Constant_Value_a[2];          /* Expression: [0.92 0]
                                        * Referenced by: '<S20>/Constant'
                                        */
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S60>/dq'
                                        */
  real_T dq_Y0_p[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S61>/dq'
                                        */
  real_T SFunction_P1_Size_aw[2];      /* Computed Parameter: SFunction_P1_Size_aw
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction_P1_p;               /* Expression: dest
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction_P2_Size_p[2];       /* Computed Parameter: SFunction_P2_Size_p
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction_P2_b;               /* Expression: priority2
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction_P3_Size_b[2];       /* Computed Parameter: SFunction_P3_Size_b
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction_P3_n;               /* Expression: st
                                        * Referenced by: '<S90>/S-Function'
                                        */
  real_T SFunction1_Value;             /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  real_T SFunction_X0;                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
                                        */
  real_T Memory_1_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Memory_2_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Memory_3_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T OpMonitor_P1_Size[2];         /* Computed Parameter: OpMonitor_P1_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P1;                 /* Expression: compute_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P2_Size[2];         /* Computed Parameter: OpMonitor_P2_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P2;                 /* Expression: real_step
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P3_Size[2];         /* Computed Parameter: OpMonitor_P3_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P3;                 /* Expression: idle_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P4_Size[2];         /* Computed Parameter: OpMonitor_P4_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P4;                 /* Expression: nb_overruns
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P5_Size[2];         /* Computed Parameter: OpMonitor_P5_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P5;                 /* Expression: user_time
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P6_Size[2];         /* Computed Parameter: OpMonitor_P6_Size
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T OpMonitor_P6[32];             /* Computed Parameter: OpMonitor_P6
                                        * Referenced by: '<S2>/OpMonitor'
                                        */
  real_T SFunction_P1_Size_m[2];       /* Computed Parameter: SFunction_P1_Size_m
                                        * Referenced by: '<S89>/S-Function'
                                        */
  real_T SFunction_P1_m;               /* Expression: Acqu_group
                                        * Referenced by: '<S89>/S-Function'
                                        */
  real_T Constant1_Value;              /* Expression: 8.55
                                        * Referenced by: '<S2>/Constant1'
                                        */
  real_T Constant2_Value;              /* Expression: 37.4*14
                                        * Referenced by: '<S2>/Constant2'
                                        */
  real_T Constant3_Value;              /* Expression: 0.06
                                        * Referenced by: '<S2>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: 60*14
                                        * Referenced by: '<S2>/Constant4'
                                        */
  real_T SFunction_P1_Size_o[2];       /* Computed Parameter: SFunction_P1_Size_o
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T SFunction_P1_g;               /* Expression: src
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T SFunction_P2_Size_m[2];       /* Computed Parameter: SFunction_P2_Size_m
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T SFunction_P2_k;               /* Expression: Data_width
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T SFunction_P3_Size_a[2];       /* Computed Parameter: SFunction_P3_Size_a
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T SFunction_P3_l;               /* Expression: st
                                        * Referenced by: '<S91>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_a; /* Expression: sps.Finit
                                        * Referenced by: '<S31>/Unit Delay'
                                        */
  real_T avoiddivisionbyzero_UpperSat; /* Expression: 70
                                        * Referenced by: '<S20>/avoid division by zero'
                                        */
  real_T avoiddivisionbyzero_LowerSat; /* Expression: 40
                                        * Referenced by: '<S20>/avoid division by zero'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 1/4
                                        * Referenced by: '<S20>/Gain'
                                        */
  real_T SFunction_P1_Size_f[2];       /* Computed Parameter: SFunction_P1_Size_f
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P1_k;               /* Expression: MaxDelay
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P2_Size_j[2];       /* Computed Parameter: SFunction_P2_Size_j
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P2_g;               /* Expression: Ts
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P3_Size_h5[2];      /* Computed Parameter: SFunction_P3_Size_h5
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P3_c;               /* Expression: InitialValue
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P4_Size_o[2];       /* Computed Parameter: SFunction_P4_Size_o
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P4_j;               /* Expression: DFT
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S31>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S31>/Discrete-Time Integrator'
                                        */
  real_T Constant4_Value_i;            /* Expression: 2*pi
                                        * Referenced by: '<S31>/Constant4'
                                        */
  real_T FirstcycleofsimulationId092Iq0_;/* Expression: 0
                                          * Referenced by: '<S20>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T FirstcycleofsimulationId092Iq_o;/* Expression: 1
                                          * Referenced by: '<S20>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T Duk_Gain;                     /* Expression: sps.D
                                        * Referenced by: '<S52>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition[2]; /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S52>/Delay_x1'
                                        */
  real_T C11_Gain;                     /* Expression: sps.C11
                                        * Referenced by: '<S55>/C11'
                                        */
  real_T Delay_x2_InitialCondition[2]; /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S52>/Delay_x2'
                                        */
  real_T C12_Gain;                     /* Expression: sps.C12
                                        * Referenced by: '<S55>/C12'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay2'
                                        */
  real_T Iq_ref_Value;                 /* Expression: 0
                                        * Referenced by: '<S4>/Iq_ref'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S15>/Integrator'
                                        */
  real_T Integrator_IC;                /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S15>/Integrator'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: 2
                                        * Referenced by: '<S16>/Gain1'
                                        */
  real_T Integ4_gainval_o;             /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S23>/Integ4'
                                        */
  real_T Integ4_IC_a;                  /* Expression: 0
                                        * Referenced by: '<S23>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_j;/* Expression: 1e6
                                          * Referenced by: '<S23>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_l;/* Expression: eps
                                          * Referenced by: '<S23>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_fw[2];      /* Computed Parameter: SFunction_P1_Size_fw
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P1_a2;              /* Expression: MaxDelay
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P2_Size_a[2];       /* Computed Parameter: SFunction_P2_Size_a
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P2_e;               /* Expression: Ts
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P3_Size_l[2];       /* Computed Parameter: SFunction_P3_Size_l
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P3_i;               /* Expression: InitialValue
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P4_Size_n[2];       /* Computed Parameter: SFunction_P4_Size_n
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P4_b;               /* Expression: DFT
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_cj;/* Expression: 0
                                        * Referenced by: '<S24>/Unit Delay'
                                        */
  real_T Constant_Value_p;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S23>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_p;/* Expression: sps.Vinit
                                        * Referenced by: '<S23>/Unit Delay1'
                                        */
  real_T Gain3_Gain_i;                 /* Expression: 2
                                        * Referenced by: '<S16>/Gain3'
                                        */
  real_T Integ4_gainval_f;             /* Computed Parameter: Integ4_gainval_f
                                        * Referenced by: '<S26>/Integ4'
                                        */
  real_T Integ4_IC_l;                  /* Expression: 0
                                        * Referenced by: '<S26>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_c;/* Expression: 1e6
                                          * Referenced by: '<S26>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_b;/* Expression: eps
                                          * Referenced by: '<S26>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_l[2];       /* Computed Parameter: SFunction_P1_Size_l
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P1_f;               /* Expression: MaxDelay
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P2_Size_f[2];       /* Computed Parameter: SFunction_P2_Size_f
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P2_d;               /* Expression: Ts
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P3_Size_f[2];       /* Computed Parameter: SFunction_P3_Size_f
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P3_h;               /* Expression: InitialValue
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P4_Size_c[2];       /* Computed Parameter: SFunction_P4_Size_c
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T SFunction_P4_g;               /* Expression: DFT
                                        * Referenced by: '<S28>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_l; /* Expression: 0
                                        * Referenced by: '<S27>/Unit Delay'
                                        */
  real_T Constant_Value_m;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S26>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_k;/* Expression: sps.Vinit
                                        * Referenced by: '<S26>/Unit Delay1'
                                        */
  real_T RadDeg_Gain_f;                /* Expression: 180/pi
                                        * Referenced by: '<S16>/Rad->Deg.'
                                        */
  real_T torad_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S11>/to-rad'
                                        */
  real_T Rff_Gain;                     /* Expression: RLff(1)
                                        * Referenced by: '<S8>/Rff '
                                        */
  real_T Lff_Gain;                     /* Expression: RLff(2)
                                        * Referenced by: '<S8>/Lff  '
                                        */
  real_T Rff_Gain_a;                   /* Expression: RLff(1)
                                        * Referenced by: '<S8>/Rff'
                                        */
  real_T Lff_Gain_c;                   /* Expression: RLff(2)
                                        * Referenced by: '<S8>/Lff'
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: 1.5
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T Saturation_LowerSat_k;        /* Expression: -1.5
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T Iph_1_Value;                  /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S9>/Iph_1'
                                        */
  real_T Iph_2_Value;                  /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S9>/Iph_2'
                                        */
  real_T MPPT_On_Value;                /* Expression: 1
                                        * Referenced by: '<S2>/MPPT_On'
                                        */
  real_T Integ4_gainval_j;             /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S29>/Integ4'
                                        */
  real_T Integ4_IC_i;                  /* Expression: 0
                                        * Referenced by: '<S29>/Integ4'
                                        */
  real_T K1_Value;                     /* Expression: sps.Delay
                                        * Referenced by: '<S29>/K1'
                                        */
  real_T SFunction_P1_Size_p[2];       /* Computed Parameter: SFunction_P1_Size_p
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P1_km;              /* Expression: MaxDelay
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P2_Size_ar[2];      /* Computed Parameter: SFunction_P2_Size_ar
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P2_ge;              /* Expression: Ts
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P3_Size_e[2];       /* Computed Parameter: SFunction_P3_Size_e
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P3_p;               /* Expression: InitialValue
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P4_Size_c0[2];      /* Computed Parameter: SFunction_P4_Size_c0
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T SFunction_P4_f;               /* Expression: DFT
                                        * Referenced by: '<S30>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S29>/Unit Delay'
                                        */
  real_T K2_Value;                     /* Expression: sps.Freq
                                        * Referenced by: '<S29>/K2'
                                        */
  real_T UnitDelay1_InitialCondition_kn;/* Expression: sps.Vinit
                                         * Referenced by: '<S29>/Unit Delay1'
                                         */
  real_T Integ4_gainval_jy;            /* Computed Parameter: Integ4_gainval_jy
                                        * Referenced by: '<S45>/Integ4'
                                        */
  real_T Integ4_IC_h;                  /* Expression: 0
                                        * Referenced by: '<S45>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_e;/* Expression: 1e6
                                          * Referenced by: '<S45>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_a;/* Expression: eps
                                          * Referenced by: '<S45>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_ol[2];      /* Computed Parameter: SFunction_P1_Size_ol
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P1_h;               /* Expression: MaxDelay
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P2_Size_h[2];       /* Computed Parameter: SFunction_P2_Size_h
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P2_h;               /* Expression: Ts
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P3_Size_j[2];       /* Computed Parameter: SFunction_P3_Size_j
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P3_nb;              /* Expression: InitialValue
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P4_Size_d[2];       /* Computed Parameter: SFunction_P4_Size_d
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T SFunction_P4_a;               /* Expression: DFT
                                        * Referenced by: '<S47>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S46>/Unit Delay'
                                        */
  real_T Constant_Value_f;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S45>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_j;/* Expression: sps.Vinit
                                        * Referenced by: '<S45>/Unit Delay1'
                                        */
  real_T DiscreteDerivative_DenCoef[2];/* Expression: [ TcD  Ts-TcD ]
                                        * Referenced by: '<S33>/Discrete Derivative '
                                        */
  real_T DiscreteDerivative_InitialState;/* Expression: 0
                                          * Referenced by: '<S33>/Discrete Derivative '
                                          */
  real_T DiscreteTimeIntegrator_gainva_j;/* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                                          * Referenced by: '<S33>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperSat;/* Expression: Par_Limits(1)
                                          * Referenced by: '<S33>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerSat;/* Expression: Par_Limits(2)
                                          * Referenced by: '<S33>/Discrete-Time Integrator'
                                          */
  real_T Saturation1_UpperSat;         /* Expression: Par_Limits(1)
                                        * Referenced by: '<S33>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: Par_Limits(2)
                                        * Referenced by: '<S33>/Saturation1'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1/2/pi
                                        * Referenced by: '<S31>/Gain10'
                                        */
  real_T RateLimiter_RisingLim;        /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S31>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S31>/Rate Limiter'
                                        */
  real_T RateLimiter_IC;               /* Expression: sps.Finit
                                        * Referenced by: '<S31>/Rate Limiter'
                                        */
  real_T Delay_x1_InitialCondition_i;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S48>/Delay_x1'
                                        */
  real_T A11_Gain;                     /* Expression: sps.A11
                                        * Referenced by: '<S49>/A11'
                                        */
  real_T Delay_x2_InitialCondition_g;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S48>/Delay_x2'
                                        */
  real_T A12_Gain;                     /* Expression: sps.A12
                                        * Referenced by: '<S49>/A12'
                                        */
  real_T A21_Gain;                     /* Expression: sps.A21
                                        * Referenced by: '<S49>/A21'
                                        */
  real_T A22_Gain;                     /* Expression: sps.A22
                                        * Referenced by: '<S49>/A22'
                                        */
  real_T B11_Gain;                     /* Expression: sps.B11
                                        * Referenced by: '<S50>/B11'
                                        */
  real_T B21_Gain;                     /* Expression: sps.B21
                                        * Referenced by: '<S50>/B21'
                                        */
  real_T Duk_Gain_p;                   /* Expression: sps.D
                                        * Referenced by: '<S48>/D*u(k)'
                                        */
  real_T C11_Gain_l;                   /* Expression: sps.C11
                                        * Referenced by: '<S51>/C11'
                                        */
  real_T C12_Gain_k;                   /* Expression: sps.C12
                                        * Referenced by: '<S51>/C12'
                                        */
  real_T A11_Gain_g;                   /* Expression: sps.A11
                                        * Referenced by: '<S53>/A11'
                                        */
  real_T A12_Gain_l;                   /* Expression: sps.A12
                                        * Referenced by: '<S53>/A12'
                                        */
  real_T A21_Gain_k;                   /* Expression: sps.A21
                                        * Referenced by: '<S53>/A21'
                                        */
  real_T A22_Gain_p;                   /* Expression: sps.A22
                                        * Referenced by: '<S53>/A22'
                                        */
  real_T B11_Gain_n;                   /* Expression: sps.B11
                                        * Referenced by: '<S54>/B11'
                                        */
  real_T B21_Gain_m;                   /* Expression: sps.B21
                                        * Referenced by: '<S54>/B21'
                                        */
  real_T Constant1_Value_p;            /* Expression: 0
                                        * Referenced by: '<S20>/Constant1'
                                        */
  real_T Constant3_Value_c;            /* Expression: sps.Delay
                                        * Referenced by: '<S84>/Constant3'
                                        */
  real_T Constant1_Value_h;            /* Expression: sps.Period
                                        * Referenced by: '<S84>/Constant1'
                                        */
  real_T ib1_Gain;                     /* Expression: sps.Freq
                                        * Referenced by: '<S84>/1\ib1'
                                        */
  real_T LookupTable_XData[3];         /* Expression: [0 .5 1]
                                        * Referenced by: '<S84>/Lookup Table'
                                        */
  real_T LookupTable_YData[3];         /* Expression: [0 2 0]
                                        * Referenced by: '<S84>/Lookup Table'
                                        */
  real_T Constant2_Value_a;            /* Expression: 1
                                        * Referenced by: '<S84>/Constant2'
                                        */
  real_T Gain1_Gain_du;                /* Expression: 0.5
                                        * Referenced by: '<S62>/Gain1'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0.1684
                                        * Referenced by: '<S4>/Unit Delay'
                                        */
  real_T Gain_Gain_cv;                 /* Expression: -1
                                        * Referenced by: '<S67>/Gain'
                                        */
  real_T Constant2_Value_b;            /* Expression: 0
                                        * Referenced by: '<S13>/Constant2'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay3'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 1
                                        * Referenced by: '<S13>/Gain1'
                                        */
  real_T UnitDelay1_InitialCondition_n;/* Expression: 0
                                        * Referenced by: '<S4>/Unit Delay1'
                                        */
  real_T Integrator_gainval_i;         /* Computed Parameter: Integrator_gainval_i
                                        * Referenced by: '<S86>/Integrator'
                                        */
  real_T Integrator_IC_i;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S86>/Integrator'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_my_pv_system_3_ss_controller_T {
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
    SimStruct childSFunctions[11];
    SimStruct *childSFunctionPtrs[11];
    struct _ssBlkInfo2 blkInfo2[11];
    struct _ssSFcnModelMethods2 methods2[11];
    struct _ssSFcnModelMethods3 methods3[11];
    struct _ssStatesInfo2 statesInfo2[11];
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
      real_T const *UPtrs0[6];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn2;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[1];
      struct _ssPortOutputs outputPortInfo[4];
      uint_T attribs[6];
      mxArray *params[6];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[4];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn4;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
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
    } Sfcn9;

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
    } Sfcn10;
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
extern P_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_P;

/* Block signals (auto storage) */
extern B_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_B;

/* Block states (auto storage) */
extern DW_my_pv_system_3_ss_controller_T my_pv_system_3_ss_controller_DW;

/* External data declarations for dependent source files */
extern const real_T my_pv_system_3_ss_controller_RGND;/* real_T ground */

/*====================*
 * External functions *
 *====================*/
extern my_pv_system_3_ss_controller_rtModel *my_pv_system_3_ss_controller(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_my_pv_system_3_ss_controller_T *const
  my_pv_system_3_ss_controller_M;

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
 * '<Root>' : 'my_pv_system_3_ss_controller'
 * '<S1>'   : 'my_pv_system_3_ss_controller/OpCCode_do_not_touch'
 * '<S2>'   : 'my_pv_system_3_ss_controller/SS_controller'
 * '<S3>'   : 'my_pv_system_3_ss_controller/powergui'
 * '<S4>'   : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control'
 * '<S5>'   : 'my_pv_system_3_ss_controller/SS_controller/MATLAB Function'
 * '<S6>'   : 'my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem'
 * '<S7>'   : 'my_pv_system_3_ss_controller/SS_controller/zzzOpComm'
 * '<S8>'   : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/Current Regulator'
 * '<S9>'   : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/MPPT  Parameters'
 * '<S10>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/MPPT Controller using Perturbe  & Observe technique  '
 * '<S11>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements'
 * '<S12>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator'
 * '<S13>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/U_ref Generation '
 * '<S14>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/VDC Regulator'
 * '<S15>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/Current Regulator/PI'
 * '<S16>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)'
 * '<S17>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean'
 * '<S18>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL'
 * '<S19>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter'
 * '<S20>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform'
 * '<S21>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S22>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S23>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S24>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S25>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S26>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S27>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S28>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S29>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean/Model'
 * '<S30>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Mean/Model/Discrete Variable Time Delay'
 * '<S31>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model'
 * '<S32>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control'
 * '<S33>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Discrete'
 * '<S34>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)'
 * '<S35>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter'
 * '<S36>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)'
 * '<S37>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S38>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S39>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S40>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S41>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S42>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S43>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S44>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S45>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model'
 * '<S46>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Correction subsystem'
 * '<S47>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Discrete Variable Time Delay'
 * '<S48>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model'
 * '<S49>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/A*k(k-1)'
 * '<S50>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S51>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/C*x(k)'
 * '<S52>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model'
 * '<S53>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/A*k(k-1)'
 * '<S54>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S55>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Second-Order Filter/Model/C*x(k)'
 * '<S56>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0'
 * '<S57>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Discrete Variable Time Delay'
 * '<S58>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S59>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S60>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S61>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S62>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Cr_MinMax'
 * '<S63>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type'
 * '<S64>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal'
 * '<S65>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling'
 * '<S66>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/Full Bridge Bipolar'
 * '<S67>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/Full Bridge Unipolar'
 * '<S68>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Modulator type/One Three Phase Bridge'
 * '<S69>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal/External'
 * '<S70>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Reference signal/Internal'
 * '<S71>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical'
 * '<S72>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Natural'
 * '<S73>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical'
 * '<S74>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical'
 * '<S75>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural'
 * '<S76>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical'
 * '<S77>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical/Sample & Hold'
 * '<S78>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Natural/Sync_NaturalSampling'
 * '<S79>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling'
 * '<S80>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling/Sample & Hold'
 * '<S81>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical/Unsync_AsymmetricalSampling'
 * '<S82>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling'
 * '<S83>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator'
 * '<S84>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator/Model'
 * '<S85>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical/Unsync_SymmetricalSampling'
 * '<S86>'  : 'my_pv_system_3_ss_controller/SS_controller/Inverter Control/VDC Regulator/PI'
 * '<S87>'  : 'my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem2'
 * '<S88>'  : 'my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem4'
 * '<S89>'  : 'my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem2/Send2'
 * '<S90>'  : 'my_pv_system_3_ss_controller/SS_controller/rtlab_send_subsystem/Subsystem4/Send4'
 * '<S91>'  : 'my_pv_system_3_ss_controller/SS_controller/zzzOpComm/Receive_1'
 * '<S92>'  : 'my_pv_system_3_ss_controller/SS_controller/zzzOpComm/busStruct'
 * '<S93>'  : 'my_pv_system_3_ss_controller/SS_controller/zzzOpComm/busStruct/Sub1'
 * '<S94>'  : 'my_pv_system_3_ss_controller/SS_controller/zzzOpComm/busStruct/Sub2'
 */
#endif                                 /* RTW_HEADER_my_pv_system_3_ss_controller_h_ */
