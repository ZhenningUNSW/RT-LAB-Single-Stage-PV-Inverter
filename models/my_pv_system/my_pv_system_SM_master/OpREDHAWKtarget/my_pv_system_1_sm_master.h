/*
 * my_pv_system_1_sm_master.h
 *
 * Code generation for model "my_pv_system_1_sm_master".
 *
 * Model version              : 1.189
 * Simulink Coder version : 8.7 (R2014b) 08-Sep-2014
 * C source code generated on : Tue Sep 06 11:04:11 2016
 *
 * Target selection: rtlab_rtmodel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#ifndef RTW_HEADER_my_pv_system_1_sm_master_h_
#define RTW_HEADER_my_pv_system_1_sm_master_h_
#include <stddef.h>
#include <string.h>
#include <math.h>
#ifndef my_pv_system_1_sm_master_COMMON_INCLUDES_
# define my_pv_system_1_sm_master_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* my_pv_system_1_sm_master_COMMON_INCLUDES_ */

#include "my_pv_system_1_sm_master_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"

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
# define rtmIsContinuousTask(rtm, tid) ((tid) == 0)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) ((rtmIsMajorTimeStep((rtm)) && (rtm)->Timing.sampleHits[(rtm)->Timing.sampleTimeTaskIDPtr[sti]]))
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
#define my_pv_system_1_sm_master_rtModel RT_MODEL_my_pv_system_1_sm_master_T

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
  real_T SFunction;                    /* '<S1>/S-Function' */
  real_T Sum;                          /* '<S1>/Sum' */
  real_T IrradianceWm2;                /* '<S2>/Memory' */
  real_T I_PV;                         /* '<S2>/Memory' */
  real_T V_PV;                         /* '<S2>/Memory' */
  real_T I_Diode;                      /* '<S2>/Memory' */
  real_T Temperature;                  /* '<S2>/Memory' */
  real_T Vph;                          /* '<S2>/Memory' */
  real_T Vinv;                         /* '<S2>/Memory1' */
  real_T Vdc;                          /* '<S2>/Memory1' */
  real_T Igrid;                        /* '<S2>/Memory1' */
  real_T Vgrid;                        /* '<S2>/Memory1' */
  real_T Ig;                           /* '<S2>/Memory1' */
  real_T OpMonitor_o1;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o2;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o3;                 /* '<S2>/OpMonitor' */
  real_T OpMonitor_o4;                 /* '<S2>/OpMonitor' */
  real_T UnitDelay[4];                 /* '<S21>/Unit Delay' */
  real_T SFunction_g[6];               /* '<S45>/S-Function' */
  real_T AC;                           /* '<S38>/AC' */
  real_T StateSpace_o1[12];            /* '<S48>/State-Space' */
  real_T StateSpace_o2[4];             /* '<S48>/State-Space' */
  real_T donotdeletethisgain;          /* '<S5>/do not delete this gain' */
  real_T donotdeletethisgain_e;        /* '<S10>/do not delete this gain' */
  real_T donotdeletethisgain_g;        /* '<S11>/do not delete this gain' */
  real_T Ron[4];                       /* '<S21>/1//Ron' */
  real_T DataTypeConversion[4];        /* '<S21>/Data Type Conversion' */
  real_T Switch[4];                    /* '<S21>/Switch' */
  real_T Saturation[4];                /* '<S21>/Saturation' */
  real_T Switch_e[4];                  /* '<S23>/Switch' */
  real_T SFunction_k[2];               /* '<S25>/S-Function' */
  real_T donotdeletethisgain_a;        /* '<S28>/do not delete this gain' */
  real_T donotdeletethisgain_d;        /* '<S29>/do not delete this gain' */
  real_T donotdeletethisgain_db;       /* '<S30>/do not delete this gain' */
  real_T donotdeletethisgain_i;        /* '<S17>/do not delete this gain' */
  real_T donotdeletethisgain_de;       /* '<S18>/do not delete this gain' */
  real_T u;                            /* '<S22>/2' */
  real_T DiscreteTimeIntegrator[4];    /* '<S22>/Discrete-Time Integrator' */
  real_T u9Tf[4];                      /* '<S24>/-0.9//Tf' */
  real_T Add[4];                       /* '<S24>/Add' */
  real_T Saturation1[4];               /* '<S24>/Saturation1' */
  real_T Add1[4];                      /* '<S24>/Add1' */
  real_T uTt[4];                       /* '<S24>/0.1//Tt' */
  real_T Saturation2[4];               /* '<S24>/Saturation2' */
  real_T Add2[4];                      /* '<S24>/Add2' */
  real_T UnitDelay_l[4];               /* '<S22>/Unit Delay' */
  real_T Switch_j[4];                  /* '<S22>/Switch' */
  real_T Product[4];                   /* '<S22>/Product' */
} B_my_pv_system_1_sm_master_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE[4];          /* '<S21>/Unit Delay' */
  real_T StateSpace_DSTATE[10];        /* '<S48>/State-Space' */
  real_T DiscreteTimeIntegrator_DSTATE[4];/* '<S22>/Discrete-Time Integrator' */
  real_T UnitDelay_DSTATE_e[4];        /* '<S22>/Unit Delay' */
  real_T SFunction_PreviousInput;      /* '<S1>/S-Function' */
  real_T Memory_1_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_2_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_3_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_4_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_5_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_6_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory1_1_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_2_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_3_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_4_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_5_PreviousInput;      /* '<S2>/Memory1' */
  void *OpMonitor_PWORK;               /* '<S2>/OpMonitor' */
  void *SFunction_PWORK;               /* '<S45>/S-Function' */
  struct {
    void *AS;
    void *BS;
    void *CS;
    void *DS;
    void *DX_COL;
    void *BD_COL;
    void *TMP1;
    void *TMP2;
    void *XTMP;
    void *SWITCH_STATUS;
    void *SWITCH_STATUS_INIT;
    void *SW_CHG;
    void *CHOPPER;
    void *G_STATE;
    void *XKM12;
    void *XKP12;
    void *XLAST;
    void *ULAST;
    void *IDX_SW_CHG;
    void *Y_SWITCH;
    void *SWITCH_TYPES;
    void *IDX_OUT_SW;
  } StateSpace_PWORK;                  /* '<S48>/State-Space' */

  int_T SFunction_IWORK[5];            /* '<S42>/S-Function' */
  int_T SFunction_IWORK_g[5];          /* '<S43>/S-Function' */
  int_T StateSpace_IWORK[5];           /* '<S48>/State-Space' */
  int_T SFunction_IWORK_b[5];          /* '<S44>/S-Function' */
  int8_T DiscreteTimeIntegrator_PrevRese[4];/* '<S22>/Discrete-Time Integrator' */
  boolean_T Tail_MODE;                 /* '<S21>/Tail' */
} DW_my_pv_system_1_sm_master_T;

/* Backward compatible GRT Identifiers */
#define rtB                            my_pv_system_1_sm_master_B
#define BlockIO                        B_my_pv_system_1_sm_master_T
#define rtP                            my_pv_system_1_sm_master_P
#define Parameters                     P_my_pv_system_1_sm_master_T
#define rtDWork                        my_pv_system_1_sm_master_DW
#define D_Work                         DW_my_pv_system_1_sm_master_T

/* Parameters (auto storage) */
struct P_my_pv_system_1_sm_master_T_ {
  real_T Vs14400V_Amplitude;           /* Mask Parameter: Vs14400V_Amplitude
                                        * Referenced by: '<S38>/AC'
                                        */
  real_T Vs14400V_Frequency;           /* Mask Parameter: Vs14400V_Frequency
                                        * Referenced by: '<S38>/AC'
                                        */
  real_T Vs14400V_Phase;               /* Mask Parameter: Vs14400V_Phase
                                        * Referenced by: '<S38>/AC'
                                        */
  real_T Tail_Tf;                      /* Mask Parameter: Tail_Tf
                                        * Referenced by:
                                        *   '<S22>/2'
                                        *   '<S24>/Constant2'
                                        *   '<S24>/-0.9//Tf'
                                        */
  real_T Tail_Tt;                      /* Mask Parameter: Tail_Tt
                                        * Referenced by:
                                        *   '<S22>/2'
                                        *   '<S24>/Constant2'
                                        *   '<S24>/0.1//Tt'
                                        */
  real_T itail_Y0;                     /* Expression: 0
                                        * Referenced by: '<S22>/itail'
                                        */
  real_T _Value;                       /* Expression: 1
                                        * Referenced by: '<S22>/1'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S22>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S22>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value;               /* Expression: 0.9
                                        * Referenced by: '<S24>/Constant'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 0.9
                                        * Referenced by: '<S24>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S24>/Saturation1'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: 0.1
                                        * Referenced by: '<S24>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S24>/Saturation2'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S22>/Unit Delay'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S22>/Switch'
                                        */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: dest
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: priority2
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S44>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: st
                                        * Referenced by: '<S44>/S-Function'
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
  real_T Memory_4_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Memory_5_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Memory_6_X0;                  /* Expression: 0
                                        * Referenced by: '<S2>/Memory'
                                        */
  real_T Memory1_1_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T Memory1_2_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T Memory1_3_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T Memory1_4_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T Memory1_5_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory1'
                                        */
  real_T SFunction_P1_Size_n[2];       /* Computed Parameter: SFunction_P1_Size_n
                                        * Referenced by: '<S42>/S-Function'
                                        */
  real_T SFunction_P1_h;               /* Expression: Acqu_group
                                        * Referenced by: '<S42>/S-Function'
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
  real_T SFunction_P1_Size_j[2];       /* Computed Parameter: SFunction_P1_Size_j
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T SFunction_P1_a;               /* Expression: Acqu_group
                                        * Referenced by: '<S43>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_p; /* Expression: 0
                                        * Referenced by: '<S21>/Unit Delay'
                                        */
  real_T SFunction_P1_Size_p[2];       /* Computed Parameter: SFunction_P1_Size_p
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T SFunction_P1_h2;              /* Expression: src
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T SFunction_P2_Size_p[2];       /* Computed Parameter: SFunction_P2_Size_p
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T SFunction_P2_f;               /* Expression: Data_width
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T SFunction_P3_Size_m[2];       /* Computed Parameter: SFunction_P3_Size_m
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T SFunction_P3_l;               /* Expression: st
                                        * Referenced by: '<S45>/S-Function'
                                        */
  real_T AC_Bias;                      /* Expression: 0
                                        * Referenced by: '<S38>/AC'
                                        */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S5>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_g;   /* Expression: 1
                                        * Referenced by: '<S10>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_h;   /* Expression: 1
                                        * Referenced by: '<S11>/do not delete this gain'
                                        */
  real_T u_Value;                      /* Expression: 0
                                        * Referenced by: '<S21>/0 4'
                                        */
  real_T Ron_Gain;                     /* Expression: 1./Ron
                                        * Referenced by: '<S21>/1//Ron'
                                        */
  real_T Switch_Threshold_l;           /* Expression: 0.5
                                        * Referenced by: '<S21>/Switch'
                                        */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S21>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0
                                        * Referenced by: '<S21>/Saturation'
                                        */
  real_T VfDevicesClampingDiodes_Value[4];/* Expression: Vf_SwitchOn
                                           * Referenced by: '<S23>/Vf Devices & Clamping Diodes'
                                           */
  real_T VfDiodes_Value[4];            /* Expression: Vf_DiodeOn
                                        * Referenced by: '<S23>/Vf Diodes'
                                        */
  real_T SFunction_P1_Size_k[2];       /* Computed Parameter: SFunction_P1_Size_k
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P1_j;               /* Expression: Data_width
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P2_Size_h[2];       /* Computed Parameter: SFunction_P2_Size_h
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T SFunction_P2_c[2];            /* Expression: InitialConditions
                                        * Referenced by: '<S25>/S-Function'
                                        */
  real_T donotdeletethisgain_Gain_j;   /* Expression: 1
                                        * Referenced by: '<S28>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_n;   /* Expression: 1
                                        * Referenced by: '<S29>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_jy;  /* Expression: 1
                                        * Referenced by: '<S30>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_c;   /* Expression: 1
                                        * Referenced by: '<S17>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_k;   /* Expression: 1
                                        * Referenced by: '<S18>/do not delete this gain'
                                        */
  boolean_T _Value_g;                  /* Expression: Tf_sps>0 | Tt_sps>0
                                        * Referenced by: '<S21>/2'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_my_pv_system_1_sm_master_T {
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
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[6];
    SimStruct *childSFunctionPtrs[6];
    struct _ssBlkInfo2 blkInfo2[6];
    struct _ssSFcnModelMethods2 methods2[6];
    struct _ssSFcnModelMethods3 methods3[6];
    struct _ssStatesInfo2 statesInfo2[6];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[11];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn0;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[11];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn1;

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
    } Sfcn2;

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
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn4;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[2];
      mxArray *params[2];
    } Sfcn5;
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
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
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
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_P;

/* Block signals (auto storage) */
extern B_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_B;

/* Block states (auto storage) */
extern DW_my_pv_system_1_sm_master_T my_pv_system_1_sm_master_DW;

/* External data declarations for dependent source files */
extern const real_T my_pv_system_1_sm_master_RGND;/* real_T ground */

/*====================*
 * External functions *
 *====================*/
extern my_pv_system_1_sm_master_rtModel *my_pv_system_1_sm_master(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_my_pv_system_1_sm_master_T *const my_pv_system_1_sm_master_M;

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
 * '<Root>' : 'my_pv_system_1_sm_master'
 * '<S1>'   : 'my_pv_system_1_sm_master/OpCCode_do_not_touch'
 * '<S2>'   : 'my_pv_system_1_sm_master/SM_master'
 * '<S3>'   : 'my_pv_system_1_sm_master/powergui'
 * '<S4>'   : 'my_pv_system_1_sm_master/SM_master/75-kVA Transformer (120V//120V//14.4 kV)'
 * '<S5>'   : 'my_pv_system_1_sm_master/SM_master/CM2'
 * '<S6>'   : 'my_pv_system_1_sm_master/SM_master/GRID'
 * '<S7>'   : 'my_pv_system_1_sm_master/SM_master/H-Bridge '
 * '<S8>'   : 'my_pv_system_1_sm_master/SM_master/OpComm'
 * '<S9>'   : 'my_pv_system_1_sm_master/SM_master/Subsystem'
 * '<S10>'  : 'my_pv_system_1_sm_master/SM_master/Vdc '
 * '<S11>'  : 'my_pv_system_1_sm_master/SM_master/Vdc 1'
 * '<S12>'  : 'my_pv_system_1_sm_master/SM_master/Vs 14,400 V'
 * '<S13>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem'
 * '<S14>'  : 'my_pv_system_1_sm_master/SM_master/zzzOpComm1'
 * '<S15>'  : 'my_pv_system_1_sm_master/SM_master/75-kVA Transformer (120V//120V//14.4 kV)/Model'
 * '<S16>'  : 'my_pv_system_1_sm_master/SM_master/CM2/Model'
 * '<S17>'  : 'my_pv_system_1_sm_master/SM_master/GRID/CM2'
 * '<S18>'  : 'my_pv_system_1_sm_master/SM_master/GRID/Vm2'
 * '<S19>'  : 'my_pv_system_1_sm_master/SM_master/GRID/CM2/Model'
 * '<S20>'  : 'my_pv_system_1_sm_master/SM_master/GRID/Vm2/Model'
 * '<S21>'  : 'my_pv_system_1_sm_master/SM_master/H-Bridge /Model'
 * '<S22>'  : 'my_pv_system_1_sm_master/SM_master/H-Bridge /Model/Tail'
 * '<S23>'  : 'my_pv_system_1_sm_master/SM_master/H-Bridge /Model/Vf 1'
 * '<S24>'  : 'my_pv_system_1_sm_master/SM_master/H-Bridge /Model/Tail/y=f(t)'
 * '<S25>'  : 'my_pv_system_1_sm_master/SM_master/OpComm/Receive'
 * '<S26>'  : 'my_pv_system_1_sm_master/SM_master/OpComm/busStruct'
 * '<S27>'  : 'my_pv_system_1_sm_master/SM_master/OpComm/busStruct/Sub1'
 * '<S28>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Current Measurement'
 * '<S29>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Voltage Measurement'
 * '<S30>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Voltage Measurement1'
 * '<S31>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/io '
 * '<S32>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/iph '
 * '<S33>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Current Measurement/Model'
 * '<S34>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Voltage Measurement/Model'
 * '<S35>'  : 'my_pv_system_1_sm_master/SM_master/Subsystem/Voltage Measurement1/Model'
 * '<S36>'  : 'my_pv_system_1_sm_master/SM_master/Vdc /Model'
 * '<S37>'  : 'my_pv_system_1_sm_master/SM_master/Vdc 1/Model'
 * '<S38>'  : 'my_pv_system_1_sm_master/SM_master/Vs 14,400 V/Model'
 * '<S39>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem1'
 * '<S40>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem2'
 * '<S41>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem5'
 * '<S42>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem1/Send1'
 * '<S43>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem2/Send2'
 * '<S44>'  : 'my_pv_system_1_sm_master/SM_master/rtlab_send_subsystem/Subsystem5/Send5'
 * '<S45>'  : 'my_pv_system_1_sm_master/SM_master/zzzOpComm1/Receive_1'
 * '<S46>'  : 'my_pv_system_1_sm_master/SM_master/zzzOpComm1/busStruct'
 * '<S47>'  : 'my_pv_system_1_sm_master/SM_master/zzzOpComm1/busStruct/Sub1'
 * '<S48>'  : 'my_pv_system_1_sm_master/powergui/EquivalentModel1'
 * '<S49>'  : 'my_pv_system_1_sm_master/powergui/EquivalentModel1/Gates'
 * '<S50>'  : 'my_pv_system_1_sm_master/powergui/EquivalentModel1/Sources'
 * '<S51>'  : 'my_pv_system_1_sm_master/powergui/EquivalentModel1/Status'
 * '<S52>'  : 'my_pv_system_1_sm_master/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_my_pv_system_1_sm_master_h_ */
