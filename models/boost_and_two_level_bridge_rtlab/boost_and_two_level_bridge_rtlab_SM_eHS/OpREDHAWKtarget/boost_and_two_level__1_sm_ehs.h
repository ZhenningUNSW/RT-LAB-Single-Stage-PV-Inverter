/*
 * boost_and_two_level__1_sm_ehs.h
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
#ifndef RTW_HEADER_boost_and_two_level__1_sm_ehs_h_
#define RTW_HEADER_boost_and_two_level__1_sm_ehs_h_
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#ifndef boost_and_two_level__1_sm_ehs_COMMON_INCLUDES_
# define boost_and_two_level__1_sm_ehs_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* boost_and_two_level__1_sm_ehs_COMMON_INCLUDES_ */

#include "boost_and_two_level__1_sm_ehs_types.h"

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
# define rtmGetIntgData(rtm)           ((rtm)->ModelData.intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->ModelData.intgData = (val))
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
# define rtmGetOdeF(rtm)               ((rtm)->ModelData.odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->ModelData.odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->ModelData.odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->ModelData.odeY = (val))
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
#define boost_and_two_level__1_sm_ehs_rtModel RT_MODEL_boost_and_two_level__1_sm_ehs_T

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

/* Block signals for system '<S6>/RMS ' */
typedef struct {
  creal_T RealImagtoComplex;           /* '<S20>/Real-Imag to Complex' */
  real_T integrator;                   /* '<S24>/integrator' */
  real_T TransportDelay;               /* '<S24>/Transport Delay' */
  real_T Clock;                        /* '<S24>/Clock' */
  real_T Memory;                       /* '<S24>/Memory' */
  real_T Switch;                       /* '<S24>/Switch' */
  real_T integrator_a;                 /* '<S23>/integrator' */
  real_T TransportDelay_p;             /* '<S23>/Transport Delay' */
  real_T Clock_h;                      /* '<S23>/Clock' */
  real_T Memory_d;                     /* '<S23>/Memory' */
  real_T Switch_e;                     /* '<S23>/Switch' */
  real_T ComplextoMagnitudeAngle_o1;   /* '<S20>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2;   /* '<S20>/Complex to Magnitude-Angle' */
  real_T sinwt;                        /* '<S20>/sin(wt)' */
  real_T Product;                      /* '<S20>/Product' */
  real_T coswt;                        /* '<S20>/cos(wt)' */
  real_T Product1;                     /* '<S20>/Product1' */
  real_T RadDeg;                       /* '<S20>/Rad->Deg.' */
  real_T Gain;                         /* '<S18>/Gain' */
  real_T Sum;                          /* '<S24>/Sum' */
  real_T Gain_l;                       /* '<S24>/Gain' */
  real_T Sum_c;                        /* '<S23>/Sum' */
  real_T Gain_n;                       /* '<S23>/Gain' */
  boolean_T RelationalOperator;        /* '<S24>/Relational Operator' */
  boolean_T RelationalOperator_a;      /* '<S23>/Relational Operator' */
} B_RMS_boost_and_two_level__1__T;

/* Block states (auto storage) for system '<S6>/RMS ' */
typedef struct {
  real_T Memory_PreviousInput;         /* '<S24>/Memory' */
  real_T Memory_PreviousInput_m;       /* '<S23>/Memory' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[16384];
  } TransportDelay_RWORK;              /* '<S24>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[16384];
  } TransportDelay_RWORK_a;            /* '<S23>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S24>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_n;            /* '<S23>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S24>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_b;            /* '<S23>/Transport Delay' */

  boolean_T RMS_MODE;                  /* '<S6>/RMS ' */
} DW_RMS_boost_and_two_level__1_T;

/* Continuous states for system '<S6>/RMS ' */
typedef struct {
  real_T integrator_CSTATE_a;          /* '<S24>/integrator' */
  real_T integrator_CSTATE_o;          /* '<S23>/integrator' */
} X_RMS_boost_and_two_level__1__T;

/* State derivatives for system '<S6>/RMS ' */
typedef struct {
  real_T integrator_CSTATE_a;          /* '<S24>/integrator' */
  real_T integrator_CSTATE_o;          /* '<S23>/integrator' */
} XDot_RMS_boost_and_two_level__T;

/* State Disabled for system '<S6>/RMS ' */
typedef struct {
  boolean_T integrator_CSTATE_a;       /* '<S24>/integrator' */
  boolean_T integrator_CSTATE_o;       /* '<S23>/integrator' */
} XDis_RMS_boost_and_two_level__T;

/* Block signals for system '<S6>/TrueRMS ' */
typedef struct {
  real_T Clock;                        /* '<S26>/Clock' */
  real_T integrator;                   /* '<S26>/integrator' */
  real_T TransportDelay;               /* '<S26>/Transport Delay' */
  real_T Memory;                       /* '<S26>/Memory' */
  real_T Switch;                       /* '<S26>/Switch' */
  real_T Product;                      /* '<S19>/Product' */
  real_T Saturationtoavoidnegativesqrt;/* '<S19>/Saturation to avoid negative sqrt' */
  real_T Sqrt;                         /* '<S19>/Sqrt' */
  real_T Sum;                          /* '<S26>/Sum' */
  real_T Gain;                         /* '<S26>/Gain' */
  boolean_T RelationalOperator;        /* '<S26>/Relational Operator' */
} B_TrueRMS_boost_and_two_level_T;

/* Block states (auto storage) for system '<S6>/TrueRMS ' */
typedef struct {
  real_T Memory_PreviousInput;         /* '<S26>/Memory' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[16384];
  } TransportDelay_RWORK;              /* '<S26>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S26>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S26>/Transport Delay' */

  boolean_T TrueRMS_MODE;              /* '<S6>/TrueRMS ' */
} DW_TrueRMS_boost_and_two_leve_T;

/* Continuous states for system '<S6>/TrueRMS ' */
typedef struct {
  real_T integrator_CSTATE_e;          /* '<S26>/integrator' */
} X_TrueRMS_boost_and_two_level_T;

/* State derivatives for system '<S6>/TrueRMS ' */
typedef struct {
  real_T integrator_CSTATE_e;          /* '<S26>/integrator' */
} XDot_TrueRMS_boost_and_two_le_T;

/* State Disabled for system '<S6>/TrueRMS ' */
typedef struct {
  boolean_T integrator_CSTATE_e;       /* '<S26>/integrator' */
} XDis_TrueRMS_boost_and_two_le_T;

/* Block signals (auto storage) */
typedef struct {
  creal_T RealImagtoComplex;           /* '<S72>/Real-Imag to Complex' */
  creal_T MagnitudeAngletoComplex;     /* '<S67>/Magnitude-Angle to Complex' */
  creal_T RealImagtoComplex_j;         /* '<S69>/Real-Imag to Complex' */
  creal_T RealImagtoComplex_f;         /* '<S92>/Real-Imag to Complex' */
  real_T SFunction;                    /* '<S1>/S-Function' */
  real_T Sum;                          /* '<S1>/Sum' */
  real_T I_LEAK;                       /* '<S2>/Memory1' */
  real_T Gate;                         /* '<S2>/Memory1' */
  real_T ext_freq;                     /* '<S2>/Memory1' */
  real_T ext_duty;                     /* '<S2>/Memory1' */
  real_T SFunction_g[165];             /* '<S15>/S-Function' */
  real_T Relay1;                       /* '<S10>/Relay1' */
  real_T DataTypeConversion8;          /* '<S10>/Data Type Conversion8' */
  real_T eHS_rst_loadin;               /* '<S10>/eHS_rst_loadin' */
  real_T Automated_Solver_Mat_Initialisa;/* '<S10>/Automated_Solver_Mat_Initialisation_1' */
  real_T SineWaveFunction;             /* '<S12>/Sine Wave Function' */
  real_T Memory1;                      /* '<S12>/Memory1' */
  real_T Memory;                       /* '<S12>/Memory' */
  real_T Inputs_eHS1_Send;             /* '<S10>/Inputs_eHS1_Send' */
  real_T Outputs_eHS1_Recv_o2;         /* '<S10>/Outputs_eHS1_Recv' */
  real_T Computation_time;             /* '<S2>/OpMonitor' */
  real_T Real_step_size;               /* '<S2>/OpMonitor' */
  real_T Idle_time;                    /* '<S2>/OpMonitor' */
  real_T Num_overruns;                 /* '<S2>/OpMonitor' */
  real_T V_INV;                        /* '<S2>/Memory2' */
  real_T I_GRID;                       /* '<S2>/Memory2' */
  real_T V_GRID;                       /* '<S2>/Memory2' */
  real_T I_PV;                         /* '<S2>/Memory' */
  real_T V_PV;                         /* '<S2>/Memory' */
  real_T P_PV;                         /* '<S2>/Memory' */
  real_T Irradiancewm2;                /* '<S2>/Memory' */
  real_T Temperature;                  /* '<S2>/Memory' */
  real_T I_DIODE;                      /* '<S2>/Memory' */
  real_T V_rms;                        /* '<S2>/Memory3' */
  real_T I_rms;                        /* '<S2>/Memory3' */
  real_T Memory1_o;                    /* '<S42>/Memory1' */
  real_T DataTypeConversion1;          /* '<S46>/Data Type Conversion1' */
  real_T DataTypeConversion;           /* '<S42>/Data Type Conversion' */
  real_T Memory2;                      /* '<S42>/Memory2' */
  real_T Add;                          /* '<S42>/Add' */
  real_T Switch;                       /* '<S42>/Switch' */
  real_T Switch1;                      /* '<S42>/Switch1' */
  real_T LoadIn;                       /* '<S40>/LoadIn' */
  real_T Outputs_eHS1_Recv_o2_j;       /* '<S2>/Outputs_eHS1_Recv' */
  real_T fpga_raw[7];                  /* '<S2>/sfp2dbl' */
  real_T Constant[9];                  /* '<S41>/Constant' */
  real_T DataInSend;                   /* '<S40>/DataIn Send' */
  real_T rtlab_io_block_o1[8];         /* '<S13>/rtlab_io_block' */
  real_T rtlab_io_block_o2[8];         /* '<S13>/rtlab_io_block' */
  real_T rtlab_io_block_o3;            /* '<S13>/rtlab_io_block' */
  real_T Memory1_n;                    /* '<S13>/Memory1' */
  real_T IOTypeSel_LoadIn;             /* '<S13>/IOTypeSel_LoadIn' */
  real_T Memory_o;                     /* '<S13>/Memory' */
  real_T rtlab_io_block_o1_c[8];       /* '<S14>/rtlab_io_block' */
  real_T rtlab_io_block_o2_a[8];       /* '<S14>/rtlab_io_block' */
  real_T rtlab_io_block_o3_f;          /* '<S14>/rtlab_io_block' */
  real_T Memory1_h;                    /* '<S14>/Memory1' */
  real_T IOTypeSel_LoadIn_j;           /* '<S14>/IOTypeSel_LoadIn' */
  real_T Memory_k;                     /* '<S14>/Memory' */
  real_T integrator;                   /* '<S38>/integrator' */
  real_T TransportDelay;               /* '<S38>/Transport Delay' */
  real_T Clock;                        /* '<S38>/Clock' */
  real_T Memory_g;                     /* '<S38>/Memory' */
  real_T Switch_f;                     /* '<S38>/Switch' */
  real_T integrator_p;                 /* '<S39>/integrator' */
  real_T TransportDelay_p;             /* '<S39>/Transport Delay' */
  real_T Clock_m;                      /* '<S39>/Clock' */
  real_T Memory_e;                     /* '<S39>/Memory' */
  real_T Switch_e;                     /* '<S39>/Switch' */
  real_T P_PV_l;                       /* '<S8>/Divide' */
  real_T ConvertSinglefloatingpointFPGAt[7];/* '<S10>/Convert  Single floating-point (FPGA)  to double' */
  real_T Divide[7];                    /* '<S10>/Divide' */
  real_T OpWriteFile_o1;               /* '<S2>/OpWriteFile' */
  real_T OpWriteFile_o2;               /* '<S2>/OpWriteFile' */
  real_T Switch_h;                     /* '<S6>/Switch' */
  real_T Switch_c;                     /* '<S7>/Switch' */
  real_T DataTypeConversion1_o;        /* '<S50>/Data Type Conversion1' */
  real_T Step;                         /* '<S12>/Step' */
  real_T Product1;                     /* '<S12>/Product1' */
  real_T TmpSignalConversionAtRTEConvers[4];
  real_T RTEConversion[4];             /* '<S12>/RTE Conversion' */
  real_T UnitDelay;                    /* '<S62>/Unit Delay' */
  real_T DigitalClock;                 /* '<S140>/Digital Clock' */
  real_T Add1;                         /* '<S140>/Add1' */
  real_T MathFunction;                 /* '<S140>/Math Function' */
  real_T ib1;                          /* '<S140>/1\ib1' */
  real_T LookupTable;                  /* '<S140>/Lookup Table' */
  real_T Add3;                         /* '<S140>/Add3' */
  real_T Add3_g;                       /* '<S118>/Add3' */
  real_T Gain1;                        /* '<S118>/Gain1' */
  real_T MUL1;                         /* '<S118>/MUL1' */
  real_T Add4;                         /* '<S118>/Add4' */
  real_T DataTypeConversion_a[4];      /* '<S68>/Data Type Conversion' */
  real_T RTEConversion1[4];            /* '<S12>/RTE Conversion1' */
  real_T RTELogicalOperator1[4];       /* '<S12>/RTE Logical Operator1' */
  real_T RTESPWM_o1;                   /* '<S12>/RTE SPWM' */
  real_T RTESPWM_o2;                   /* '<S12>/RTE SPWM' */
  real_T Switch_fr[4];                 /* '<S12>/Switch' */
  real_T RTEGround[4];                 /* '<S12>/RTE Ground' */
  real_T RTE_Conversion_1_o1[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o2[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o3[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o4[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o5[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o6[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o7[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o8[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o9[4];       /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o10[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o11[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o12[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o13[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o14[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o15[4];      /* '<S49>/RTE_Conversion_1' */
  real_T RTE_Conversion_1_o16[4];      /* '<S49>/RTE_Conversion_1' */
  real_T EventGen_eHS_1;               /* '<S49>/EventGen_eHS_1' */
  real_T RateTransition5;              /* '<S12>/Rate Transition5' */
  real_T Apu;                          /* '<S67>/A->pu' */
  real_T UnitDelay_l;                  /* '<S87>/Unit Delay' */
  real_T avoiddivisionbyzero;          /* '<S76>/avoid division by zero' */
  real_T MathFunction_i;               /* '<S76>/Math Function' */
  real_T Gain;                         /* '<S76>/Gain' */
  real_T SFunction_h;                  /* '<S113>/S-Function' */
  real_T DiscreteTimeIntegrator;       /* '<S87>/Discrete-Time Integrator' */
  real_T MathFunction_m;               /* '<S87>/Math Function' */
  real_T FirstcycleofsimulationId092Iq0;/* '<S76>/First cycle of simulation Id=0.92, Iq=0' */
  real_T Switch_b[2];                  /* '<S76>/Switch' */
  real_T Duk[2];                       /* '<S108>/D*u(k)' */
  real_T x1k[2];                       /* '<S108>/Delay_x1' */
  real_T C11[2];                       /* '<S111>/C11' */
  real_T x2k[2];                       /* '<S108>/Delay_x2' */
  real_T C12[2];                       /* '<S111>/C12' */
  real_T sum2[2];                      /* '<S111>/sum2' */
  real_T yk[2];                        /* '<S108>/C*X(k)+D*u(k)' */
  real_T UnitDelay2;                   /* '<S62>/Unit Delay2' */
  real_T Sum_n[2];                     /* '<S64>/Sum' */
  real_T ProportionalGain[2];          /* '<S71>/Proportional Gain' */
  real_T Integrator[2];                /* '<S71>/Integrator' */
  real_T Sum_j[2];                     /* '<S71>/Sum' */
  real_T Saturate[2];                  /* '<S71>/Saturate' */
  real_T RateTransition4;              /* '<S12>/Rate Transition4' */
  real_T Vpu;                          /* '<S67>/V->pu' */
  real_T TrigonometricFunction;        /* '<S72>/Trigonometric Function' */
  real_T Gain1_f;                      /* '<S72>/Gain1' */
  real_T Product1_k;                   /* '<S72>/Product1' */
  real_T Integ4;                       /* '<S79>/Integ4' */
  real_T Freq;                         /* '<S79>/To avoid division  by zero' */
  real_T Numberofsamplespercycle;      /* '<S79>/Number of samples per cycle' */
  real_T RoundingFunction;             /* '<S79>/Rounding Function' */
  real_T Delay;                        /* '<S79>/Gain' */
  real_T SFunction_f;                  /* '<S81>/S-Function' */
  real_T UnitDelay_k;                  /* '<S80>/Unit Delay' */
  real_T DigitalClock_e;               /* '<S79>/Digital  Clock' */
  real_T UnitDelay1;                   /* '<S79>/Unit Delay1' */
  real_T Switch_i;                     /* '<S79>/Switch' */
  real_T TrigonometricFunction3;       /* '<S72>/Trigonometric Function3' */
  real_T Gain3;                        /* '<S72>/Gain3' */
  real_T Product2;                     /* '<S72>/Product2' */
  real_T Integ4_b;                     /* '<S82>/Integ4' */
  real_T Freq_g;                       /* '<S82>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_h;    /* '<S82>/Number of samples per cycle' */
  real_T RoundingFunction_e;           /* '<S82>/Rounding Function' */
  real_T Delay_p;                      /* '<S82>/Gain' */
  real_T SFunction_fz;                 /* '<S84>/S-Function' */
  real_T UnitDelay_h;                  /* '<S83>/Unit Delay' */
  real_T DigitalClock_i;               /* '<S82>/Digital  Clock' */
  real_T UnitDelay1_c;                 /* '<S82>/Unit Delay1' */
  real_T Switch_cy;                    /* '<S82>/Switch' */
  real_T ComplextoMagnitudeAngle_o1;   /* '<S72>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2;   /* '<S72>/Complex to Magnitude-Angle' */
  real_T RadDeg;                       /* '<S72>/Rad->Deg.' */
  real_T torad;                        /* '<S67>/to-rad' */
  real_T ComplextoRealImag_o1;         /* '<S67>/Complex to Real-Imag' */
  real_T ComplextoRealImag_o2;         /* '<S67>/Complex to Real-Imag' */
  real_T Rff;                          /* '<S64>/Rff ' */
  real_T Lff;                          /* '<S64>/Lff  ' */
  real_T Feedforward;                  /* '<S64>/Add1' */
  real_T Rff_d;                        /* '<S64>/Rff' */
  real_T Lff_o;                        /* '<S64>/Lff' */
  real_T Add3_i;                       /* '<S64>/Add3' */
  real_T Add2[2];                      /* '<S64>/Add2' */
  real_T IntegralGain[2];              /* '<S71>/Integral Gain' */
  real_T Saturation[2];                /* '<S64>/Saturation' */
  real_T RateTransition7;              /* '<S12>/Rate Transition7' */
  real_T RateTransition8;              /* '<S12>/Rate Transition8' */
  real_T DigitalClock_j;               /* '<S85>/Digital  Clock' */
  real_T RateTransition6;              /* '<S12>/Rate Transition6' */
  real_T Integ4_m;                     /* '<S85>/Integ4' */
  real_T K1;                           /* '<S85>/K1' */
  real_T SFunction_c;                  /* '<S86>/S-Function' */
  real_T UnitDelay_b;                  /* '<S85>/Unit Delay' */
  real_T UnitDelay1_cm;                /* '<S85>/Unit Delay1' */
  real_T Switch_ej;                    /* '<S85>/Switch' */
  real_T TrigonometricFunction2;       /* '<S87>/Trigonometric Function2' */
  real_T Product1_g;                   /* '<S87>/Product1' */
  real_T Integ4_e;                     /* '<S101>/Integ4' */
  real_T Freq_a;                       /* '<S101>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_g;    /* '<S101>/Number of samples per cycle' */
  real_T RoundingFunction_l;           /* '<S101>/Rounding Function' */
  real_T Delay_pz;                     /* '<S101>/Gain' */
  real_T SFunction_i;                  /* '<S103>/S-Function' */
  real_T UnitDelay_lg;                 /* '<S102>/Unit Delay' */
  real_T DigitalClock_a;               /* '<S101>/Digital  Clock' */
  real_T UnitDelay1_i;                 /* '<S101>/Unit Delay1' */
  real_T Switch_cd;                    /* '<S101>/Switch' */
  real_T Divide_f;                     /* '<S87>/Divide' */
  real_T DiscreteDerivative;           /* '<S89>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_i;     /* '<S89>/Discrete-Time Integrator' */
  real_T Kp4;                          /* '<S89>/Kp4' */
  real_T Sum6;                         /* '<S89>/Sum6' */
  real_T Saturation1;                  /* '<S89>/Saturation1' */
  real_T Gain10;                       /* '<S87>/Gain10' */
  real_T RateLimiter;                  /* '<S87>/Rate Limiter' */
  real_T x1k_e;                        /* '<S104>/Delay_x1' */
  real_T A11;                          /* '<S105>/A11' */
  real_T x2k_h;                        /* '<S104>/Delay_x2' */
  real_T A12;                          /* '<S105>/A12' */
  real_T A21;                          /* '<S105>/A21' */
  real_T A22;                          /* '<S105>/A22' */
  real_T sum2_c;                       /* '<S105>/sum2' */
  real_T sum3;                         /* '<S105>/sum3' */
  real_T B11;                          /* '<S106>/B11' */
  real_T x1k1;                         /* '<S104>/A*x1(k) + B*u1(k) ' */
  real_T B21;                          /* '<S106>/B21' */
  real_T x2k1;                         /* '<S104>/A*x2(k) + B*u2(k)' */
  real_T Duk_c;                        /* '<S104>/D*u(k)' */
  real_T C11_f;                        /* '<S107>/C11' */
  real_T C12_o;                        /* '<S107>/C12' */
  real_T sum2_n;                       /* '<S107>/sum2' */
  real_T yk_c;                         /* '<S104>/C*X(k)+D*u(k)' */
  real_T A11_n[2];                     /* '<S109>/A11' */
  real_T A12_e[2];                     /* '<S109>/A12' */
  real_T A21_b[2];                     /* '<S109>/A21' */
  real_T A22_o[2];                     /* '<S109>/A22' */
  real_T sum2_a[2];                    /* '<S109>/sum2' */
  real_T sum3_f[2];                    /* '<S109>/sum3' */
  real_T B11_a[2];                     /* '<S110>/B11' */
  real_T x1k1_l[2];                    /* '<S108>/A*x1(k) + B*u1(k) ' */
  real_T B21_a[2];                     /* '<S110>/B21' */
  real_T x2k1_j[2];                    /* '<S108>/A*x2(k) + B*u2(k)' */
  real_T Constant1;                    /* '<S76>/Constant1' */
  real_T Switch_ex;                    /* '<S62>/Switch' */
  real_T Add1_j;                       /* '<S69>/Add1' */
  real_T UnitDelay3[2];                /* '<S62>/Unit Delay3' */
  real_T Gain1_p;                      /* '<S69>/Gain1' */
  real_T Product;                      /* '<S69>/Product' */
  real_T Product1_m[2];                /* '<S69>/Product1' */
  real_T ComplextoMagnitudeAngle_o1_k; /* '<S69>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_g; /* '<S69>/Complex to Magnitude-Angle' */
  real_T Add2_l;                       /* '<S69>/Add2' */
  real_T TrigonometricFunction_c;      /* '<S69>/Trigonometric Function' */
  real_T Product2_g;                   /* '<S69>/Product2' */
  real_T UnitDelay1_k;                 /* '<S62>/Unit Delay1' */
  real_T Sum_p;                        /* '<S70>/Sum' */
  real_T Rtot_pu2;                     /* '<S70>/Rtot_pu2' */
  real_T IntegralGain_m;               /* '<S142>/Integral Gain' */
  real_T Integrator_j;                 /* '<S142>/Integrator' */
  real_T ProportionalGain_e;           /* '<S142>/Proportional Gain' */
  real_T Sum_e;                        /* '<S142>/Sum' */
  real_T Saturate_d;                   /* '<S142>/Saturate' */
  real_T RateTransition3;              /* '<S12>/Rate Transition3' */
  real_T RateTransition1;              /* '<S12>/Rate Transition1' */
  real_T RateTransition2;              /* '<S12>/Rate Transition2' */
  real_T RTEPeriodMeter;               /* '<S12>/RTE Period Meter' */
  real_T Saturation_p;                 /* '<S12>/Saturation' */
  real_T Clock_o;                      /* '<S2>/Clock' */
  real_T Constant1_f;                  /* '<S2>/Constant1' */
  real_T DataTypeConversion_c;         /* '<S2>/Data Type Conversion' */
  real_T OpTrigger;                    /* '<S2>/OpTrigger' */
  real_T OpCtrl_o1;                    /* '<S2>/OpCtrl' */
  real_T OpCtrl_o2[4];                 /* '<S2>/OpCtrl' */
  real_T Iph;                          /* '<S12>/MATLAB Function' */
  real_T Io;                           /* '<S12>/MATLAB Function' */
  real_T Fcn;                          /* '<S117>/Fcn' */
  real_T Fcn1;                         /* '<S117>/Fcn1' */
  real_T Fcn_i;                        /* '<S116>/Fcn' */
  real_T Fcn1_k;                       /* '<S116>/Fcn1' */
  real_T Switch_cq[2];                 /* '<S112>/Switch' */
  real_T Sum1;                         /* '<S102>/Sum1' */
  real_T Sum5;                         /* '<S102>/Sum5' */
  real_T Product5;                     /* '<S102>/Product5' */
  real_T Gain1_e;                      /* '<S102>/Gain1' */
  real_T Sum4;                         /* '<S102>/Sum4' */
  real_T Product2_f;                   /* '<S102>/Product2' */
  real_T Product4;                     /* '<S102>/Product4' */
  real_T Sum7;                         /* '<S101>/Sum7' */
  real_T Meanvalue;                    /* '<S101>/Product' */
  real_T Sum5_h;                       /* '<S101>/Sum5' */
  real_T TrigonometricFunction_e;      /* '<S92>/Trigonometric Function' */
  real_T Gain1_fj;                     /* '<S92>/Gain1' */
  real_T Product1_a;                   /* '<S92>/Product1' */
  real_T Integ4_ee;                    /* '<S95>/Integ4' */
  real_T Freq_h;                       /* '<S95>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_o;    /* '<S95>/Number of samples per cycle' */
  real_T RoundingFunction_j;           /* '<S95>/Rounding Function' */
  real_T Delay_o;                      /* '<S95>/Gain' */
  real_T SFunction_b;                  /* '<S97>/S-Function' */
  real_T UnitDelay_e;                  /* '<S96>/Unit Delay' */
  real_T DigitalClock_g;               /* '<S95>/Digital  Clock' */
  real_T UnitDelay1_a;                 /* '<S95>/Unit Delay1' */
  real_T Switch_hk;                    /* '<S95>/Switch' */
  real_T TrigonometricFunction3_k;     /* '<S92>/Trigonometric Function3' */
  real_T Gain3_b;                      /* '<S92>/Gain3' */
  real_T Product2_b;                   /* '<S92>/Product2' */
  real_T Integ4_i;                     /* '<S98>/Integ4' */
  real_T Freq_c;                       /* '<S98>/To avoid division  by zero' */
  real_T Numberofsamplespercycle_e;    /* '<S98>/Number of samples per cycle' */
  real_T RoundingFunction_n;           /* '<S98>/Rounding Function' */
  real_T Delay_e;                      /* '<S98>/Gain' */
  real_T SFunction_j;                  /* '<S100>/S-Function' */
  real_T UnitDelay_m;                  /* '<S99>/Unit Delay' */
  real_T DigitalClock_d;               /* '<S98>/Digital  Clock' */
  real_T UnitDelay1_n;                 /* '<S98>/Unit Delay1' */
  real_T Switch_cf;                    /* '<S98>/Switch' */
  real_T ComplextoMagnitudeAngle_o1_c; /* '<S92>/Complex to Magnitude-Angle' */
  real_T ComplextoMagnitudeAngle_o2_k; /* '<S92>/Complex to Magnitude-Angle' */
  real_T RadDeg_l;                     /* '<S92>/Rad->Deg.' */
  real_T Saturation_o;                 /* '<S88>/Saturation' */
  real_T MathFunction_c;               /* '<S88>/Math Function' */
  real_T Sum1_h;                       /* '<S99>/Sum1' */
  real_T Sum5_n;                       /* '<S99>/Sum5' */
  real_T Product5_j;                   /* '<S99>/Product5' */
  real_T Gain1_b;                      /* '<S99>/Gain1' */
  real_T Sum4_j;                       /* '<S99>/Sum4' */
  real_T Product2_m;                   /* '<S99>/Product2' */
  real_T Product4_h;                   /* '<S99>/Product4' */
  real_T Sum7_g;                       /* '<S98>/Sum7' */
  real_T Meanvalue_b;                  /* '<S98>/Product' */
  real_T Sum5_nt;                      /* '<S98>/Sum5' */
  real_T Sum1_i;                       /* '<S96>/Sum1' */
  real_T Sum5_a;                       /* '<S96>/Sum5' */
  real_T Product5_i;                   /* '<S96>/Product5' */
  real_T Gain1_p5;                     /* '<S96>/Gain1' */
  real_T Sum4_n;                       /* '<S96>/Sum4' */
  real_T Product2_c;                   /* '<S96>/Product2' */
  real_T Product4_p;                   /* '<S96>/Product4' */
  real_T Sum7_j;                       /* '<S95>/Sum7' */
  real_T Meanvalue_h;                  /* '<S95>/Product' */
  real_T Sum5_d;                       /* '<S95>/Sum5' */
  real_T Gain1_a;                      /* '<S85>/Gain1' */
  real_T Gain_o;                       /* '<S85>/Gain' */
  real_T Correction;                   /* '<S85>/Sum1' */
  real_T Sum7_js;                      /* '<S85>/Sum7' */
  real_T Mean;                         /* '<S85>/Product' */
  real_T Sum5_c;                       /* '<S85>/Sum5' */
  real_T Sum1_e;                       /* '<S83>/Sum1' */
  real_T Sum5_p;                       /* '<S83>/Sum5' */
  real_T Product5_i4;                  /* '<S83>/Product5' */
  real_T Gain1_i;                      /* '<S83>/Gain1' */
  real_T Sum4_l;                       /* '<S83>/Sum4' */
  real_T Product2_gq;                  /* '<S83>/Product2' */
  real_T Product4_n;                   /* '<S83>/Product4' */
  real_T Sum7_m;                       /* '<S82>/Sum7' */
  real_T Meanvalue_a;                  /* '<S82>/Product' */
  real_T Sum5_an;                      /* '<S82>/Sum5' */
  real_T Sum1_g;                       /* '<S80>/Sum1' */
  real_T Sum5_f;                       /* '<S80>/Sum5' */
  real_T Product5_l;                   /* '<S80>/Product5' */
  real_T Gain1_pu;                     /* '<S80>/Gain1' */
  real_T Sum4_f;                       /* '<S80>/Sum4' */
  real_T Product2_f1;                  /* '<S80>/Product2' */
  real_T Product4_k;                   /* '<S80>/Product4' */
  real_T Sum7_p;                       /* '<S79>/Sum7' */
  real_T Meanvalue_o;                  /* '<S79>/Product' */
  real_T Sum5_n4;                      /* '<S79>/Sum5' */
  real_T TmpSignalConversionAtSFunctionI[4];/* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T D;                            /* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Saturation_f[16];             /* '<S45>/Saturation' */
  real_T Gain1_l[16];                  /* '<S45>/Gain1' */
  real_T Saturation1_p[16];            /* '<S45>/Saturation1' */
  real_T Gain2[16];                    /* '<S45>/Gain2' */
  real_T Saturation_g[16];             /* '<S44>/Saturation' */
  real_T Gain1_j[16];                  /* '<S44>/Gain1' */
  real_T Saturation1_c[16];            /* '<S44>/Saturation1' */
  real_T Gain2_d[16];                  /* '<S44>/Gain2' */
  real_T Saturation_h[16];             /* '<S43>/Saturation' */
  real_T Gain1_g[16];                  /* '<S43>/Gain1' */
  real_T Saturation1_cu[16];           /* '<S43>/Saturation1' */
  real_T Gain2_a[16];                  /* '<S43>/Gain2' */
  real_T Sum_ns;                       /* '<S39>/Sum' */
  real_T Gain_k;                       /* '<S39>/Gain' */
  real_T Sum_h;                        /* '<S38>/Sum' */
  real_T Gain_a;                       /* '<S38>/Gain' */
  int32_T DataTypeConversion2[16];     /* '<S45>/Data Type Conversion2' */
  int32_T DataTypeConversion2_c[16];   /* '<S44>/Data Type Conversion2' */
  int32_T DataTypeConversion2_n[16];   /* '<S43>/Data Type Conversion2' */
  uint32_T DataTypeConversion1_h;      /* '<S10>/Data Type Conversion1' */
  uint32_T sat_scn;                    /* '<S10>/sat_scn' */
  uint32_T shift_2bits;                /* '<S10>/shift_2bits' */
  uint32_T Add_l;                      /* '<S10>/Add' */
  uint32_T Uk1;                        /* '<S48>/Delay Input1' */
  uint32_T ConvertdoubletoSinglefloatingpo[3];/* '<S10>/Convert double to  Single floating-point (FPGA)' */
  uint32_T Outputs_eHS1_Recv_o1[8];    /* '<S10>/Outputs_eHS1_Recv' */
  uint32_T ConvertdoubletoSinglefloating_o[16];/* '<S40>/Convert double to  Single floating-point (FPGA)' */
  uint32_T ConvertdoubletoSinglefloatin_ou[16];/* '<S40>/Convert double to  Single floating-point (FPGA)1' */
  uint32_T Switch_n[16];               /* '<S45>/Switch' */
  uint32_T BitwiseOperator[16];        /* '<S45>/Bitwise Operator' */
  uint32_T DataTypeConversion5[16];    /* '<S45>/Data Type Conversion5' */
  uint32_T Switch_g[16];               /* '<S43>/Switch' */
  uint32_T BitwiseOperator_p[16];      /* '<S43>/Bitwise Operator' */
  uint32_T DataTypeConversion5_m[16];  /* '<S43>/Data Type Conversion5' */
  uint32_T Switch_d[16];               /* '<S44>/Switch' */
  uint32_T BitwiseOperator_pi[16];     /* '<S44>/Bitwise Operator' */
  uint32_T DataTypeConversion5_n[16];  /* '<S44>/Data Type Conversion5' */
  uint32_T Uk1_g[80];                  /* '<S47>/Delay Input1' */
  uint32_T MultiportSwitch[17];        /* '<S40>/Multiport Switch' */
  uint32_T TmpSignalConversionAtLoadInInpo[18];
  uint32_T Outputs_eHS1_Recv_o1_h[7];  /* '<S2>/Outputs_eHS1_Recv' */
  uint32_T ConvertdoubletoSinglefloating_e[16];/* '<S40>/Convert double to  Single floating-point (FPGA)2' */
  uint32_T IOTypeSel;                  /* '<S13>/IOTypeSel' */
  uint32_T IOTypeSel_p;                /* '<S14>/IOTypeSel' */
  uint32_T Saturation_k;               /* '<S10>/Saturation' */
  uint32_T load_config1[10];           /* '<S10>/load_config1' */
  uint32_T Uk1_gd[10];                 /* '<S51>/Delay Input1' */
  uint32_T DataTypeConversion1_l[16];  /* '<S45>/Data Type Conversion1' */
  uint32_T DataTypeConversion1_d[16];  /* '<S44>/Data Type Conversion1' */
  uint32_T DataTypeConversion1_hm[16]; /* '<S43>/Data Type Conversion1' */
  uint8_T FixPtRelationalOperator;     /* '<S48>/FixPt Relational Operator' */
  uint8_T FixPtRelationalOperator_i[80];/* '<S47>/FixPt Relational Operator' */
  uint8_T FixPtRelationalOperator_h[10];/* '<S51>/FixPt Relational Operator' */
  uint8_T Compare;                     /* '<S114>/Compare' */
  uint8_T Compare_m;                   /* '<S115>/Compare' */
  boolean_T LogicalOperator;           /* '<S10>/Logical Operator' */
  boolean_T LogicalOperator_i;         /* '<S46>/Logical Operator' */
  boolean_T LogicalOperator1;          /* '<S42>/Logical Operator1' */
  boolean_T RelationalOperator;        /* '<S38>/Relational Operator' */
  boolean_T RelationalOperator_b;      /* '<S39>/Relational Operator' */
  boolean_T DataTypeConversion_ch;     /* '<S6>/Data Type Conversion' */
  boolean_T LogicalOperator_m;         /* '<S6>/Logical Operator' */
  boolean_T DataTypeConversion_g;      /* '<S7>/Data Type Conversion' */
  boolean_T LogicalOperator_ia;        /* '<S7>/Logical Operator' */
  boolean_T LogicalOperator_g;         /* '<S50>/Logical Operator' */
  boolean_T RelationalOperator2;       /* '<S122>/Relational Operator2' */
  boolean_T LogicalOperator_mq;        /* '<S122>/Logical Operator' */
  boolean_T LogicalOperator4[2];       /* '<S68>/Logical Operator4' */
  boolean_T RelationalOperator_m;      /* '<S79>/Relational Operator' */
  boolean_T RelationalOperator_a;      /* '<S82>/Relational Operator' */
  boolean_T RelationalOperator_h;      /* '<S85>/Relational Operator' */
  boolean_T RelationalOperator_k;      /* '<S101>/Relational Operator' */
  boolean_T Compare_p;                 /* '<S3>/Compare' */
  boolean_T RelationalOperator_mw;     /* '<S95>/Relational Operator' */
  boolean_T RelationalOperator_am;     /* '<S98>/Relational Operator' */
  B_TrueRMS_boost_and_two_level_T TrueRMS_i;/* '<S7>/TrueRMS ' */
  B_RMS_boost_and_two_level__1__T RMS_h;/* '<S7>/RMS ' */
  B_TrueRMS_boost_and_two_level_T TrueRMS;/* '<S6>/TrueRMS ' */
  B_RMS_boost_and_two_level__1__T RMS; /* '<S6>/RMS ' */
} B_boost_and_two_level__1_sm_ehs_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<S62>/Unit Delay' */
  real_T UnitDelay_DSTATE_b;           /* '<S87>/Unit Delay' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S87>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE[2];           /* '<S108>/Delay_x1' */
  real_T Delay_x2_DSTATE[2];           /* '<S108>/Delay_x2' */
  real_T UnitDelay2_DSTATE;            /* '<S62>/Unit Delay2' */
  real_T Integrator_DSTATE[2];         /* '<S71>/Integrator' */
  real_T Integ4_DSTATE;                /* '<S79>/Integ4' */
  real_T UnitDelay_DSTATE_d;           /* '<S80>/Unit Delay' */
  real_T UnitDelay1_DSTATE;            /* '<S79>/Unit Delay1' */
  real_T Integ4_DSTATE_f;              /* '<S82>/Integ4' */
  real_T UnitDelay_DSTATE_i;           /* '<S83>/Unit Delay' */
  real_T UnitDelay1_DSTATE_i;          /* '<S82>/Unit Delay1' */
  real_T Integ4_DSTATE_g;              /* '<S85>/Integ4' */
  real_T UnitDelay_DSTATE_h;           /* '<S85>/Unit Delay' */
  real_T UnitDelay1_DSTATE_im;         /* '<S85>/Unit Delay1' */
  real_T Integ4_DSTATE_k;              /* '<S101>/Integ4' */
  real_T UnitDelay_DSTATE_f;           /* '<S102>/Unit Delay' */
  real_T UnitDelay1_DSTATE_e;          /* '<S101>/Unit Delay1' */
  real_T DiscreteDerivative_states;    /* '<S89>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_DSTATE_o;/* '<S89>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE_p;            /* '<S104>/Delay_x1' */
  real_T Delay_x2_DSTATE_l;            /* '<S104>/Delay_x2' */
  real_T UnitDelay3_DSTATE[2];         /* '<S62>/Unit Delay3' */
  real_T UnitDelay1_DSTATE_n;          /* '<S62>/Unit Delay1' */
  real_T Integrator_DSTATE_h;          /* '<S142>/Integrator' */
  real_T Integ4_DSTATE_i;              /* '<S95>/Integ4' */
  real_T UnitDelay_DSTATE_c;           /* '<S96>/Unit Delay' */
  real_T UnitDelay1_DSTATE_eb;         /* '<S95>/Unit Delay1' */
  real_T Integ4_DSTATE_c;              /* '<S98>/Integ4' */
  real_T UnitDelay_DSTATE_fv;          /* '<S99>/Unit Delay' */
  real_T UnitDelay1_DSTATE_k;          /* '<S98>/Unit Delay1' */
  real_T SFunction_PreviousInput;      /* '<S1>/S-Function' */
  real_T Memory1_1_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_2_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_3_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_4_PreviousInput;      /* '<S2>/Memory1' */
  real_T Memory1_PreviousInput;        /* '<S12>/Memory1' */
  real_T Memory_PreviousInput;         /* '<S12>/Memory' */
  real_T Memory2_1_PreviousInput;      /* '<S2>/Memory2' */
  real_T Memory2_2_PreviousInput;      /* '<S2>/Memory2' */
  real_T Memory2_3_PreviousInput;      /* '<S2>/Memory2' */
  real_T Memory_1_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_2_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_3_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_4_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_5_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory_6_PreviousInput;       /* '<S2>/Memory' */
  real_T Memory3_1_PreviousInput;      /* '<S2>/Memory3' */
  real_T Memory3_2_PreviousInput;      /* '<S2>/Memory3' */
  real_T Memory1_PreviousInput_h;      /* '<S42>/Memory1' */
  real_T Memory2_PreviousInput;        /* '<S42>/Memory2' */
  real_T Memory1_PreviousInput_n;      /* '<S13>/Memory1' */
  real_T Memory_PreviousInput_l;       /* '<S13>/Memory' */
  real_T Memory1_PreviousInput_o;      /* '<S14>/Memory1' */
  real_T Memory_PreviousInput_h;       /* '<S14>/Memory' */
  real_T Memory_PreviousInput_n;       /* '<S38>/Memory' */
  real_T Memory_PreviousInput_e;       /* '<S39>/Memory' */
  real_T DiscreteDerivative_tmp;       /* '<S89>/Discrete Derivative ' */
  real_T PrevY;                        /* '<S87>/Rate Limiter' */
  real_T Vold;                         /* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Pold;                         /* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  real_T Dold;                         /* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[16384];
  } TransportDelay_RWORK;              /* '<S38>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[16384];
  } TransportDelay_RWORK_o;            /* '<S39>/Transport Delay' */

  real_T RTEConversion_RWORK[4];       /* '<S12>/RTE Conversion' */
  real_T RTEConversion1_RWORK[4];      /* '<S12>/RTE Conversion1' */
  real_T SFunction_RWORK;              /* '<S113>/S-Function' */
  real_T SFunction_RWORK_d;            /* '<S81>/S-Function' */
  real_T SFunction_RWORK_g;            /* '<S84>/S-Function' */
  real_T SFunction_RWORK_m;            /* '<S86>/S-Function' */
  real_T SFunction_RWORK_m4;           /* '<S103>/S-Function' */
  real_T SFunction_RWORK_e;            /* '<S97>/S-Function' */
  real_T SFunction_RWORK_p;            /* '<S100>/S-Function' */
  void *eHS_rst_loadin_PWORK;          /* '<S10>/eHS_rst_loadin' */
  void *Automated_Solver_Mat_Initialisa[2];/* '<S10>/Automated_Solver_Mat_Initialisation_1' */
  void *Inputs_eHS1_Send_PWORK;        /* '<S10>/Inputs_eHS1_Send' */
  void *Outputs_eHS1_Recv_PWORK;       /* '<S10>/Outputs_eHS1_Recv' */
  void *OpMonitor_PWORK;               /* '<S2>/OpMonitor' */
  void *LoadIn_PWORK;                  /* '<S40>/LoadIn' */
  void *Outputs_eHS1_Recv_PWORK_d;     /* '<S2>/Outputs_eHS1_Recv' */
  void *DataInSend_PWORK;              /* '<S40>/DataIn Send' */
  void *rtlab_io_block_PWORK;          /* '<S13>/rtlab_io_block' */
  void *IOTypeSel_LoadIn_PWORK;        /* '<S13>/IOTypeSel_LoadIn' */
  void *rtlab_io_block_PWORK_f;        /* '<S14>/rtlab_io_block' */
  void *IOTypeSel_LoadIn_PWORK_j;      /* '<S14>/IOTypeSel_LoadIn' */
  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S38>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_d;            /* '<S39>/Transport Delay' */

  void *OpWriteFile_PWORK;             /* '<S2>/OpWriteFile' */
  void *RTEConversion_PWORK;           /* '<S12>/RTE Conversion' */
  void *RTEConversion1_PWORK;          /* '<S12>/RTE Conversion1' */
  void *RTELogicalOperator1_PWORK;     /* '<S12>/RTE Logical Operator1' */
  void *RTESPWM_PWORK;                 /* '<S12>/RTE SPWM' */
  void *RTEGround_PWORK;               /* '<S12>/RTE Ground' */
  void *RTE_Conversion_1_PWORK;        /* '<S49>/RTE_Conversion_1' */
  void *EventGen_eHS_1_PWORK;          /* '<S49>/EventGen_eHS_1' */
  void *SFunction_PWORK;               /* '<S113>/S-Function' */
  struct {
    void *LoggedData;
  } PI_Ireg1_PWORK;                    /* '<S64>/PI_Ireg1' */

  void *SFunction_PWORK_o;             /* '<S81>/S-Function' */
  void *SFunction_PWORK_e;             /* '<S84>/S-Function' */
  void *SFunction_PWORK_b;             /* '<S86>/S-Function' */
  void *SFunction_PWORK_l;             /* '<S103>/S-Function' */
  void *RTEPeriodMeter_PWORK;          /* '<S12>/RTE Period Meter' */
  void *OpCtrl_PWORK;                  /* '<S2>/OpCtrl' */
  void *SFunction_PWORK_lr;            /* '<S97>/S-Function' */
  void *SFunction_PWORK_j;             /* '<S100>/S-Function' */
  uint32_T DelayInput1_DSTATE;         /* '<S48>/Delay Input1' */
  uint32_T DelayInput1_DSTATE_m[80];   /* '<S47>/Delay Input1' */
  uint32_T DelayInput1_DSTATE_k[10];   /* '<S51>/Delay Input1' */
  int_T SFunction_IWORK[5];            /* '<S57>/S-Function' */
  int_T SFunction_IWORK_b[5];          /* '<S58>/S-Function' */
  int_T SFunction_IWORK_j[5];          /* '<S59>/S-Function' */
  int_T SFunction_IWORK_o[5];          /* '<S60>/S-Function' */
  int_T SFunction_IWORK_m[5];          /* '<S61>/S-Function' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S38>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_g;            /* '<S39>/Transport Delay' */

  int_T SFunction_IWORK_n;             /* '<S113>/S-Function' */
  int_T SFunction_IWORK_bd;            /* '<S81>/S-Function' */
  int_T SFunction_IWORK_jw;            /* '<S84>/S-Function' */
  int_T SFunction_IWORK_os;            /* '<S86>/S-Function' */
  int_T SFunction_IWORK_k;             /* '<S103>/S-Function' */
  int_T OpTrigger_IWORK[5];            /* '<S2>/OpTrigger' */
  int_T OpCtrl_IWORK;                  /* '<S2>/OpCtrl' */
  int_T SFunction_IWORK_nb;            /* '<S97>/S-Function' */
  int_T SFunction_IWORK_l;             /* '<S100>/S-Function' */
  uint8_T Integ4_SYSTEM_ENABLE;        /* '<S79>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_n;      /* '<S82>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_i;      /* '<S85>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_e;      /* '<S101>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_nz;     /* '<S95>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_d;      /* '<S98>/Integ4' */
  boolean_T Relay1_Mode;               /* '<S10>/Relay1' */
  boolean_T Vold_not_empty;            /* '<S62>/MPPT Controller using Perturbe  & Observe technique  ' */
  boolean_T AutomaticGainControl_MODE; /* '<S87>/Automatic Gain Control' */
  DW_TrueRMS_boost_and_two_leve_T TrueRMS_i;/* '<S7>/TrueRMS ' */
  DW_RMS_boost_and_two_level__1_T RMS_h;/* '<S7>/RMS ' */
  DW_TrueRMS_boost_and_two_leve_T TrueRMS;/* '<S6>/TrueRMS ' */
  DW_RMS_boost_and_two_level__1_T RMS; /* '<S6>/RMS ' */
} DW_boost_and_two_level__1_sm_ehs_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T integrator_CSTATE;            /* '<S38>/integrator' */
  real_T integrator_CSTATE_b;          /* '<S39>/integrator' */
  X_TrueRMS_boost_and_two_level_T TrueRMS_i;/* '<S6>/TrueRMS ' */
  X_RMS_boost_and_two_level__1__T RMS_h;/* '<S6>/RMS ' */
  X_TrueRMS_boost_and_two_level_T TrueRMS;/* '<S6>/TrueRMS ' */
  X_RMS_boost_and_two_level__1__T RMS; /* '<S6>/RMS ' */
} X_boost_and_two_level__1_sm_ehs_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T integrator_CSTATE;            /* '<S38>/integrator' */
  real_T integrator_CSTATE_b;          /* '<S39>/integrator' */
  XDot_TrueRMS_boost_and_two_le_T TrueRMS_i;/* '<S6>/TrueRMS ' */
  XDot_RMS_boost_and_two_level__T RMS_h;/* '<S6>/RMS ' */
  XDot_TrueRMS_boost_and_two_le_T TrueRMS;/* '<S6>/TrueRMS ' */
  XDot_RMS_boost_and_two_level__T RMS; /* '<S6>/RMS ' */
} XDot_boost_and_two_level__1_sm_ehs_T;

/* State disabled  */
typedef struct {
  boolean_T integrator_CSTATE;         /* '<S38>/integrator' */
  boolean_T integrator_CSTATE_b;       /* '<S39>/integrator' */
  XDis_TrueRMS_boost_and_two_le_T TrueRMS_i;/* '<S6>/TrueRMS ' */
  XDis_RMS_boost_and_two_level__T RMS_h;/* '<S6>/RMS ' */
  XDis_TrueRMS_boost_and_two_le_T TrueRMS;/* '<S6>/TrueRMS ' */
  XDis_RMS_boost_and_two_level__T RMS; /* '<S6>/RMS ' */
} XDis_boost_and_two_level__1_sm_ehs_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Backward compatible GRT Identifiers */
#define rtB                            boost_and_two_level__1_sm_ehs_B
#define BlockIO                        B_boost_and_two_level__1_sm_ehs_T
#define rtX                            boost_and_two_level__1_sm_ehs_X
#define ContinuousStates               X_boost_and_two_level__1_sm_ehs_T
#define rtXdot                         boost_and_two_level__1_sm_ehs_XDot
#define StateDerivatives               XDot_boost_and_two_level__1_sm_ehs_T
#define tXdis                          boost_and_two_level__1_sm_ehs_XDis
#define StateDisabled                  XDis_boost_and_two_level__1_sm_ehs_T
#define rtP                            boost_and_two_level__1_sm_ehs_P
#define Parameters                     P_boost_and_two_level__1_sm_ehs_T
#define rtDWork                        boost_and_two_level__1_sm_ehs_DW
#define D_Work                         DW_boost_and_two_level__1_sm_ehs_T

/* Parameters for system: '<S6>/RMS ' */
struct P_RMS_boost_and_two_level__1__T_ {
  real_T Fourier1_Freq;                /* Mask Parameter: Fourier1_Freq
                                        * Referenced by:
                                        *   '<S20>/cos(wt)'
                                        *   '<S20>/sin(wt)'
                                        */
  real_T Gain_Gain;                    /* Expression: sps.Freq
                                        * Referenced by: '<S23>/Gain'
                                        */
  real_T Gain_Gain_l;                  /* Expression: sps.Freq
                                        * Referenced by: '<S24>/Gain'
                                        */
  real_T integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S24>/integrator'
                                        */
  real_T TransportDelay_Delay;         /* Expression: 1./sps.Freq
                                        * Referenced by: '<S24>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S24>/Transport Delay'
                                        */
  real_T K1_Value;                     /* Expression: 1./sps.Freq
                                        * Referenced by: '<S24>/K1'
                                        */
  real_T Memory_X0;                    /* Expression: sps.Vinit
                                        * Referenced by: '<S24>/Memory'
                                        */
  real_T integrator_IC_g;              /* Expression: 0
                                        * Referenced by: '<S23>/integrator'
                                        */
  real_T TransportDelay_Delay_e;       /* Expression: 1./sps.Freq
                                        * Referenced by: '<S23>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput_e;  /* Expression: 0
                                        * Referenced by: '<S23>/Transport Delay'
                                        */
  real_T K1_Value_c;                   /* Expression: 1./sps.Freq
                                        * Referenced by: '<S23>/K1'
                                        */
  real_T Memory_X0_d;                  /* Expression: sps.Vinit
                                        * Referenced by: '<S23>/Memory'
                                        */
  real_T sinwt_Amp;                    /* Expression: sps.k
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
  real_T sinwt_Bias;                   /* Expression: 0
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
  real_T sinwt_Phase;                  /* Expression: 0
                                        * Referenced by: '<S20>/sin(wt)'
                                        */
  real_T coswt_Amp;                    /* Expression: sps.k
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
  real_T coswt_Bias;                   /* Expression: 0
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
  real_T coswt_Phase;                  /* Expression: pi/2
                                        * Referenced by: '<S20>/cos(wt)'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S20>/Rad->Deg.'
                                        */
  real_T Gain_Gain_b;                  /* Expression: 1/sqrt(2)
                                        * Referenced by: '<S18>/Gain'
                                        */
};

/* Parameters for system: '<S6>/TrueRMS ' */
struct P_TrueRMS_boost_and_two_level_T_ {
  real_T Gain_Gain;                    /* Expression: sps.Freq
                                        * Referenced by: '<S26>/Gain'
                                        */
  real_T integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S26>/integrator'
                                        */
  real_T TransportDelay_Delay;         /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S26>/Transport Delay'
                                        */
  real_T K1_Value;                     /* Expression: 1./sps.Freq
                                        * Referenced by: '<S26>/K1'
                                        */
  real_T Memory_X0;                    /* Expression: sps.Vinit
                                        * Referenced by: '<S26>/Memory'
                                        */
  real_T Saturationtoavoidnegativesqrt_U;/* Expression: Inf
                                          * Referenced by: '<S19>/Saturation to avoid negative sqrt'
                                          */
  real_T Saturationtoavoidnegativesqrt_L;/* Expression: 0
                                          * Referenced by: '<S19>/Saturation to avoid negative sqrt'
                                          */
};

/* Parameters (auto storage) */
struct P_boost_and_two_level__1_sm_ehs_T_ {
  real_T Ts;                           /* Variable: Ts
                                        * Referenced by:
                                        *   '<S69>/Constant4'
                                        *   '<S79>/Gain'
                                        *   '<S82>/Gain'
                                        *   '<S101>/Gain'
                                        *   '<S95>/Gain'
                                        *   '<S98>/Gain'
                                        */
  real_T PLL_AGC;                      /* Mask Parameter: PLL_AGC
                                        * Referenced by: '<S87>/Constant1'
                                        */
  real_T AlphaBetaZerotodq0_Alignment; /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S112>/Constant'
                                        */
  real_T InverterControl_Fnom;         /* Mask Parameter: InverterControl_Fnom
                                        * Referenced by:
                                        *   '<S69>/Constant4'
                                        *   '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                        */
  real_T InverterControl_Increment_MPPT;/* Mask Parameter: InverterControl_Increment_MPPT
                                         * Referenced by: '<S65>/Iph_3'
                                         */
  real_T Discrete_Init;                /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  real_T Discrete_Kd;                  /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S89>/Discrete Derivative '
                                        */
  real_T InverterControl_Ki_Ireg;      /* Mask Parameter: InverterControl_Ki_Ireg
                                        * Referenced by: '<S71>/Integral Gain'
                                        */
  real_T InverterControl_Ki_VDCreg;    /* Mask Parameter: InverterControl_Ki_VDCreg
                                        * Referenced by: '<S142>/Integral Gain'
                                        */
  real_T Discrete_Kp;                  /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S89>/Kp4'
                                        */
  real_T InverterControl_Kp_Ireg;      /* Mask Parameter: InverterControl_Kp_Ireg
                                        * Referenced by: '<S71>/Proportional Gain'
                                        */
  real_T InverterControl_Kp_VDCreg;    /* Mask Parameter: InverterControl_Kp_VDCreg
                                        * Referenced by: '<S142>/Proportional Gain'
                                        */
  real_T PI_LowerSaturationLimit;      /* Mask Parameter: PI_LowerSaturationLimit
                                        * Referenced by: '<S71>/Saturate'
                                        */
  real_T PI_LowerSaturationLimit_j;    /* Mask Parameter: PI_LowerSaturationLimit_j
                                        * Referenced by: '<S142>/Saturate'
                                        */
  real_T PWM_Generator_MinMax[2];      /* Mask Parameter: PWM_Generator_MinMax
                                        * Referenced by: '<S68>/Constant10'
                                        */
  real_T InverterControl_Pnom;         /* Mask Parameter: InverterControl_Pnom
                                        * Referenced by: '<S67>/A->pu'
                                        */
  real_T RMS_TrueRMS;                  /* Mask Parameter: RMS_TrueRMS
                                        * Referenced by: '<S6>/Constant'
                                        */
  real_T RMS1_TrueRMS;                 /* Mask Parameter: RMS1_TrueRMS
                                        * Referenced by: '<S7>/Constant'
                                        */
  real_T PI_UpperSaturationLimit;      /* Mask Parameter: PI_UpperSaturationLimit
                                        * Referenced by: '<S71>/Saturate'
                                        */
  real_T PI_UpperSaturationLimit_a;    /* Mask Parameter: PI_UpperSaturationLimit_a
                                        * Referenced by: '<S142>/Saturate'
                                        */
  real_T InverterControl_Vdc_ref_Init; /* Mask Parameter: InverterControl_Vdc_ref_Init
                                        * Referenced by:
                                        *   '<S62>/Vnom_dc1'
                                        *   '<S65>/Iph_'
                                        */
  real_T InverterControl_Vnom_dc;      /* Mask Parameter: InverterControl_Vnom_dc
                                        * Referenced by: '<S70>/Rtot_pu2'
                                        */
  real_T InverterControl_Vnom_prim;    /* Mask Parameter: InverterControl_Vnom_prim
                                        * Referenced by:
                                        *   '<S67>/A->pu'
                                        *   '<S67>/V->pu'
                                        *   '<S69>/Constant3'
                                        */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S114>/Constant'
                                        */
  real_T CompareToConstant1_const;     /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S115>/Constant'
                                        */
  real_T CompareToConstant_const_d;    /* Mask Parameter: CompareToConstant_const_d
                                        * Referenced by: '<S3>/Constant'
                                        */
  real_T Counter_max_count;            /* Mask Parameter: Counter_max_count
                                        * Referenced by: '<S42>/Switch'
                                        */
  uint32_T BitwiseOperator_BitMask;    /* Mask Parameter: BitwiseOperator_BitMask
                                        * Referenced by: '<S45>/Bitwise Operator'
                                        */
  uint32_T BitwiseOperator_BitMask_j;  /* Mask Parameter: BitwiseOperator_BitMask_j
                                        * Referenced by: '<S43>/Bitwise Operator'
                                        */
  uint32_T BitwiseOperator_BitMask_h;  /* Mask Parameter: BitwiseOperator_BitMask_h
                                        * Referenced by: '<S44>/Bitwise Operator'
                                        */
  uint32_T DetectChange_vinit;         /* Mask Parameter: DetectChange_vinit
                                        * Referenced by: '<S48>/Delay Input1'
                                        */
  uint32_T DetectChange_vinit_i;       /* Mask Parameter: DetectChange_vinit_i
                                        * Referenced by: '<S47>/Delay Input1'
                                        */
  uint32_T DetectChange_vinit_g;       /* Mask Parameter: DetectChange_vinit_g
                                        * Referenced by: '<S51>/Delay Input1'
                                        */
  real_T Gain_Gain;                    /* Expression: sps.Freq
                                        * Referenced by: '<S38>/Gain'
                                        */
  real_T Gain_Gain_l;                  /* Expression: sps.Freq
                                        * Referenced by: '<S39>/Gain'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S42>/Constant1'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S42>/Constant2'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S42>/Constant'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S43>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S43>/Saturation1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 2^Q
                                        * Referenced by: '<S43>/Gain2'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S43>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S43>/Saturation'
                                        */
  real_T Gain1_Gain;                   /* Expression: 2^Q
                                        * Referenced by: '<S43>/Gain1'
                                        */
  real_T Saturation1_UpperSat_k;       /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S44>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_n;       /* Expression: 0
                                        * Referenced by: '<S44>/Saturation1'
                                        */
  real_T Gain2_Gain_d;                 /* Expression: 2^Q
                                        * Referenced by: '<S44>/Gain2'
                                        */
  real_T Saturation_UpperSat_i;        /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S44>/Saturation'
                                        */
  real_T Saturation_LowerSat_h;        /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S44>/Saturation'
                                        */
  real_T Gain1_Gain_j;                 /* Expression: 2^Q
                                        * Referenced by: '<S44>/Gain1'
                                        */
  real_T Saturation1_UpperSat_f;       /* Expression: 2^(N-Q)-2^(-Q)
                                        * Referenced by: '<S45>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_a;       /* Expression: 0
                                        * Referenced by: '<S45>/Saturation1'
                                        */
  real_T Gain2_Gain_p;                 /* Expression: 2^Q
                                        * Referenced by: '<S45>/Gain2'
                                        */
  real_T Saturation_UpperSat_d;        /* Expression: 2^(N-Q-1)-2^(-Q)
                                        * Referenced by: '<S45>/Saturation'
                                        */
  real_T Saturation_LowerSat_g;        /* Expression: -2^(N-Q-1)
                                        * Referenced by: '<S45>/Saturation'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 2^Q
                                        * Referenced by: '<S45>/Gain1'
                                        */
  real_T Gain1_Gain_jj;                /* Expression: 0.5
                                        * Referenced by: '<S80>/Gain1'
                                        */
  real_T Gain1_Gain_jb;                /* Expression: 0.5
                                        * Referenced by: '<S83>/Gain1'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: sps.K2
                                        * Referenced by: '<S85>/Gain1'
                                        */
  real_T Gain_Gain_j;                  /* Expression: sps.K1
                                        * Referenced by: '<S85>/Gain'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 0.5
                                        * Referenced by: '<S96>/Gain1'
                                        */
  real_T Gain1_Gain_e3;                /* Expression: 0.5
                                        * Referenced by: '<S99>/Gain1'
                                        */
  real_T Gain_Y0;                      /* Expression: [1]
                                        * Referenced by: '<S88>/Gain'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: 2
                                        * Referenced by: '<S92>/Gain1'
                                        */
  real_T Integ4_gainval;               /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S95>/Integ4'
                                        */
  real_T Integ4_IC;                    /* Expression: 0
                                        * Referenced by: '<S95>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSat;/* Expression: 1e6
                                         * Referenced by: '<S95>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_LowerSat;/* Expression: eps
                                         * Referenced by: '<S95>/To avoid division  by zero'
                                         */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: MaxDelay
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: Ts
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: InitialValue
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P4_Size[2];         /* Computed Parameter: SFunction_P4_Size
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T SFunction_P4;                 /* Expression: DFT
                                        * Referenced by: '<S97>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S96>/Unit Delay'
                                        */
  real_T Constant_Value_l;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S95>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: sps.Vinit
                                        * Referenced by: '<S95>/Unit Delay1'
                                        */
  real_T Gain3_Gain;                   /* Expression: 2
                                        * Referenced by: '<S92>/Gain3'
                                        */
  real_T Integ4_gainval_c;             /* Computed Parameter: Integ4_gainval_c
                                        * Referenced by: '<S98>/Integ4'
                                        */
  real_T Integ4_IC_b;                  /* Expression: 0
                                        * Referenced by: '<S98>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_i;/* Expression: 1e6
                                          * Referenced by: '<S98>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_h;/* Expression: eps
                                          * Referenced by: '<S98>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_h[2];       /* Computed Parameter: SFunction_P1_Size_h
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P1_i;               /* Expression: MaxDelay
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P2_Size_p[2];       /* Computed Parameter: SFunction_P2_Size_p
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P2_i;               /* Expression: Ts
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P3_Size_p[2];       /* Computed Parameter: SFunction_P3_Size_p
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P3_b;               /* Expression: InitialValue
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P4_Size_p[2];       /* Computed Parameter: SFunction_P4_Size_p
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T SFunction_P4_d;               /* Expression: DFT
                                        * Referenced by: '<S100>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_k; /* Expression: 0
                                        * Referenced by: '<S99>/Unit Delay'
                                        */
  real_T Constant_Value_o;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S98>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_n;/* Expression: sps.Vinit
                                        * Referenced by: '<S98>/Unit Delay1'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S92>/Rad->Deg.'
                                        */
  real_T Saturation_UpperSat_b;        /* Expression: inf
                                        * Referenced by: '<S88>/Saturation'
                                        */
  real_T Saturation_LowerSat_n;        /* Expression: eps
                                        * Referenced by: '<S88>/Saturation'
                                        */
  real_T Gain1_Gain_n;                 /* Expression: 0.5
                                        * Referenced by: '<S102>/Gain1'
                                        */
  real_T Constant_Value_k[2];          /* Expression: [0.92 0]
                                        * Referenced by: '<S76>/Constant'
                                        */
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S116>/dq'
                                        */
  real_T dq_Y0_b[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S117>/dq'
                                        */
  real_T SFunction1_Value;             /* Expression: 0
                                        * Referenced by: '<S1>/S-Function1'
                                        */
  real_T SFunction_X0;                 /* Expression: 0
                                        * Referenced by: '<S1>/S-Function'
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
  real_T SFunction_P1_Size_m[2];       /* Computed Parameter: SFunction_P1_Size_m
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P1_k;               /* Expression: Acqu_group
                                        * Referenced by: '<S57>/S-Function'
                                        */
  real_T SFunction_P1_Size_n[2];       /* Computed Parameter: SFunction_P1_Size_n
                                        * Referenced by: '<S15>/S-Function'
                                        */
  real_T SFunction_P1_k2;              /* Expression: Data_width
                                        * Referenced by: '<S15>/S-Function'
                                        */
  real_T SFunction_P2_Size_n[2];       /* Computed Parameter: SFunction_P2_Size_n
                                        * Referenced by: '<S15>/S-Function'
                                        */
  real_T SFunction_P2_d[165];          /* Expression: InitialConditions
                                        * Referenced by: '<S15>/S-Function'
                                        */
  real_T Relay1_OnVal;                 /* Expression: .5
                                        * Referenced by: '<S10>/Relay1'
                                        */
  real_T Relay1_OffVal;                /* Expression: .5
                                        * Referenced by: '<S10>/Relay1'
                                        */
  real_T Relay1_YOn;                   /* Expression: 1
                                        * Referenced by: '<S10>/Relay1'
                                        */
  real_T Relay1_YOff;                  /* Expression: 0
                                        * Referenced by: '<S10>/Relay1'
                                        */
  real_T eHS_rst_loadin_P1_Size[2];    /* Computed Parameter: eHS_rst_loadin_P1_Size
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P1[6];         /* Computed Parameter: eHS_rst_loadin_P1
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P2_Size[2];    /* Computed Parameter: eHS_rst_loadin_P2_Size
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P2;            /* Expression: FcnNos
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P3_Size[2];    /* Computed Parameter: eHS_rst_loadin_P3_Size
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P3;            /* Expression: width
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P4_Size[2];    /* Computed Parameter: eHS_rst_loadin_P4_Size
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T eHS_rst_loadin_P4;            /* Expression: portType
                                        * Referenced by: '<S10>/eHS_rst_loadin'
                                        */
  real_T Automated_Solver_Mat_Initialisa[2];/* Computed Parameter: Automated_Solver_Mat_Initialisa
                                             * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                             */
  real_T Automated_Solver_Mat_Initiali_c[6];/* Computed Parameter: Automated_Solver_Mat_Initiali_c
                                             * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                             */
  real_T Automated_Solver_Mat_Initiali_o[2];/* Computed Parameter: Automated_Solver_Mat_Initiali_o
                                             * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                             */
  real_T Automated_Solver_Mat_Initiali_b;/* Expression: fpga_port_in
                                          * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                          */
  real_T Automated_Solver_Mat_Initiali_e[2];/* Computed Parameter: Automated_Solver_Mat_Initiali_e
                                             * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                             */
  real_T Automated_Solver_Mat_Initial_er[31];/* Computed Parameter: Automated_Solver_Mat_Initial_er
                                              * Referenced by: '<S10>/Automated_Solver_Mat_Initialisation_1'
                                              */
  real_T SineWaveFunction_Amp;         /* Expression: 14400*sqrt(2)
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  real_T SineWaveFunction_Bias;        /* Expression: 0
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  real_T SineWaveFunction_Freq;        /* Expression: 50*2*pi
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  real_T SineWaveFunction_Phase;       /* Expression: 0
                                        * Referenced by: '<S12>/Sine Wave Function'
                                        */
  real_T Memory1_X0;                   /* Expression: 0
                                        * Referenced by: '<S12>/Memory1'
                                        */
  real_T Memory_X0;                    /* Expression: 0
                                        * Referenced by: '<S12>/Memory'
                                        */
  real_T Inputs_eHS1_Send_P1_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P1_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P1[6];       /* Computed Parameter: Inputs_eHS1_Send_P1
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P2_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P2_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P2;          /* Expression: FcnNos
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P3_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P3_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P3;          /* Expression: width
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P4_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P4_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P4;          /* Expression: portType
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P5_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P5_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P5;          /* Expression: checkVersion
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P6_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P6_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P6;          /* Expression: expectedId
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P7_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P7_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P7;          /* Expression: expectedVersion
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P8_Size[2];  /* Computed Parameter: Inputs_eHS1_Send_P8_Size
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Inputs_eHS1_Send_P8;          /* Expression: opComp
                                        * Referenced by: '<S10>/Inputs_eHS1_Send'
                                        */
  real_T Outputs_eHS1_Recv_P1_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P1_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P1[6];      /* Computed Parameter: Outputs_eHS1_Recv_P1
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P2_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P2_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P2;         /* Expression: FcnNos
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P3_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P3_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P3;         /* Expression: width
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P4_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P4_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P4;         /* Expression: portType
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P5_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P5_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P5;         /* Expression: sampleTime
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P6_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P6_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P6;         /* Expression: checkVersion
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P7_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P7_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P7;         /* Expression: expectedId
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P8_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P8_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P8;         /* Expression: expectedVersion
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P9_Size[2]; /* Computed Parameter: Outputs_eHS1_Recv_P9_Size
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P9;         /* Expression: opComp
                                        * Referenced by: '<S10>/Outputs_eHS1_Recv'
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
  real_T SFunction_P1_Size_g[2];       /* Computed Parameter: SFunction_P1_Size_g
                                        * Referenced by: '<S58>/S-Function'
                                        */
  real_T SFunction_P1_a;               /* Expression: Acqu_group
                                        * Referenced by: '<S58>/S-Function'
                                        */
  real_T Memory2_1_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */
  real_T Memory2_2_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */
  real_T Memory2_3_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory2'
                                        */
  real_T SFunction_P1_Size_b[2];       /* Computed Parameter: SFunction_P1_Size_b
                                        * Referenced by: '<S59>/S-Function'
                                        */
  real_T SFunction_P1_e;               /* Expression: Acqu_group
                                        * Referenced by: '<S59>/S-Function'
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
  real_T SFunction_P1_Size_b0[2];      /* Computed Parameter: SFunction_P1_Size_b0
                                        * Referenced by: '<S60>/S-Function'
                                        */
  real_T SFunction_P1_c;               /* Expression: Acqu_group
                                        * Referenced by: '<S60>/S-Function'
                                        */
  real_T Memory3_1_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory3'
                                        */
  real_T Memory3_2_X0;                 /* Expression: 0
                                        * Referenced by: '<S2>/Memory3'
                                        */
  real_T SFunction_P1_Size_k[2];       /* Computed Parameter: SFunction_P1_Size_k
                                        * Referenced by: '<S61>/S-Function'
                                        */
  real_T SFunction_P1_if;              /* Expression: Acqu_group
                                        * Referenced by: '<S61>/S-Function'
                                        */
  real_T Memory1_X0_c;                 /* Expression: 0
                                        * Referenced by: '<S42>/Memory1'
                                        */
  real_T Memory2_X0;                   /* Expression: 0
                                        * Referenced by: '<S42>/Memory2'
                                        */
  real_T LoadIn_P1_Size[2];            /* Computed Parameter: LoadIn_P1_Size
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P1[6];                 /* Computed Parameter: LoadIn_P1
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P2_Size[2];            /* Computed Parameter: LoadIn_P2_Size
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P2;                    /* Expression: FcnNos
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P3_Size[2];            /* Computed Parameter: LoadIn_P3_Size
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P3;                    /* Expression: width
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P4_Size[2];            /* Computed Parameter: LoadIn_P4_Size
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T LoadIn_P4;                    /* Expression: portType
                                        * Referenced by: '<S40>/LoadIn'
                                        */
  real_T Outputs_eHS1_Recv_P1_Size_e[2];/* Computed Parameter: Outputs_eHS1_Recv_P1_Size_e
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P1_p[6];    /* Computed Parameter: Outputs_eHS1_Recv_P1_p
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P2_Size_g[2];/* Computed Parameter: Outputs_eHS1_Recv_P2_Size_g
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P2_b;       /* Expression: FcnNos
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P3_Size_a[2];/* Computed Parameter: Outputs_eHS1_Recv_P3_Size_a
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P3_a;       /* Expression: width
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P4_Size_b[2];/* Computed Parameter: Outputs_eHS1_Recv_P4_Size_b
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P4_l;       /* Expression: portType
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P5_Size_m[2];/* Computed Parameter: Outputs_eHS1_Recv_P5_Size_m
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P5_j;       /* Expression: sampleTime
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P6_Size_b[2];/* Computed Parameter: Outputs_eHS1_Recv_P6_Size_b
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P6_o;       /* Expression: checkVersion
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P7_Size_j[2];/* Computed Parameter: Outputs_eHS1_Recv_P7_Size_j
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P7_i;       /* Expression: expectedId
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P8_Size_p[2];/* Computed Parameter: Outputs_eHS1_Recv_P8_Size_p
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P8_f;       /* Expression: expectedVersion
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Outputs_eHS1_Recv_P9_Size_m[2];/* Computed Parameter: Outputs_eHS1_Recv_P9_Size_m
                                         * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                         */
  real_T Outputs_eHS1_Recv_P9_i;       /* Expression: opComp
                                        * Referenced by: '<S2>/Outputs_eHS1_Recv'
                                        */
  real_T Constant_Value_g[9];          /* Expression: ones(1,9)
                                        * Referenced by: '<S41>/Constant'
                                        */
  real_T DataInSend_P1_Size[2];        /* Computed Parameter: DataInSend_P1_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P1[6];             /* Computed Parameter: DataInSend_P1
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P2_Size[2];        /* Computed Parameter: DataInSend_P2_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P2;                /* Expression: FcnNos
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P3_Size[2];        /* Computed Parameter: DataInSend_P3_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P3;                /* Expression: width
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P4_Size[2];        /* Computed Parameter: DataInSend_P4_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P4;                /* Expression: portType
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P5_Size[2];        /* Computed Parameter: DataInSend_P5_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P5;                /* Expression: checkVersion
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P6_Size[2];        /* Computed Parameter: DataInSend_P6_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P6;                /* Expression: expectedId
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P7_Size[2];        /* Computed Parameter: DataInSend_P7_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P7;                /* Expression: expectedVersion
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P8_Size[2];        /* Computed Parameter: DataInSend_P8_Size
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T DataInSend_P8;                /* Expression: opComp
                                        * Referenced by: '<S40>/DataIn Send'
                                        */
  real_T rtlab_io_block_P1_Size[2];    /* Computed Parameter: rtlab_io_block_P1_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P1[6];         /* Computed Parameter: rtlab_io_block_P1
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P2_Size[2];    /* Computed Parameter: rtlab_io_block_P2_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P3_Size[2];    /* Computed Parameter: rtlab_io_block_P3_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P3;            /* Expression: sampleTime
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P4_Size[2];    /* Computed Parameter: rtlab_io_block_P4_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P4;            /* Expression: portNb
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P5_Size[2];    /* Computed Parameter: rtlab_io_block_P5_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P5;            /* Expression: numchan
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P6_Size[2];    /* Computed Parameter: rtlab_io_block_P6_Size
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P6;            /* Expression: maxCount
                                        * Referenced by: '<S13>/rtlab_io_block'
                                        */
  real_T IOTypeSel1_Value;             /* Expression: 0
                                        * Referenced by: '<S13>/IOTypeSel1'
                                        */
  real_T Memory1_X0_a;                 /* Expression: 1
                                        * Referenced by: '<S13>/Memory1'
                                        */
  real_T IOTypeSel_LoadIn_P1_Size[2];  /* Computed Parameter: IOTypeSel_LoadIn_P1_Size
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P1[6];       /* Computed Parameter: IOTypeSel_LoadIn_P1
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P2_Size[2];  /* Computed Parameter: IOTypeSel_LoadIn_P2_Size
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P2;          /* Expression: FcnNos
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P3_Size[2];  /* Computed Parameter: IOTypeSel_LoadIn_P3_Size
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P3;          /* Expression: width
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P4_Size[2];  /* Computed Parameter: IOTypeSel_LoadIn_P4_Size
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P4;          /* Expression: portType
                                        * Referenced by: '<S13>/IOTypeSel_LoadIn'
                                        */
  real_T Memory_X0_k;                  /* Expression: 1
                                        * Referenced by: '<S13>/Memory'
                                        */
  real_T rtlab_io_block_P1_Size_b[2];  /* Computed Parameter: rtlab_io_block_P1_Size_b
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P1_m[6];       /* Computed Parameter: rtlab_io_block_P1_m
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P2_Size_k[2];  /* Computed Parameter: rtlab_io_block_P2_Size_k
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P3_Size_i[2];  /* Computed Parameter: rtlab_io_block_P3_Size_i
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P3_c;          /* Expression: sampleTime
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P4_Size_g[2];  /* Computed Parameter: rtlab_io_block_P4_Size_g
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P4_b;          /* Expression: portNb
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P5_Size_j[2];  /* Computed Parameter: rtlab_io_block_P5_Size_j
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P5_p;          /* Expression: numchan
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P6_Size_k[2];  /* Computed Parameter: rtlab_io_block_P6_Size_k
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T rtlab_io_block_P6_m;          /* Expression: maxCount
                                        * Referenced by: '<S14>/rtlab_io_block'
                                        */
  real_T IOTypeSel1_Value_d;           /* Expression: 0
                                        * Referenced by: '<S14>/IOTypeSel1'
                                        */
  real_T Memory1_X0_o;                 /* Expression: 1
                                        * Referenced by: '<S14>/Memory1'
                                        */
  real_T IOTypeSel_LoadIn_P1_Size_e[2];/* Computed Parameter: IOTypeSel_LoadIn_P1_Size_e
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P1_l[6];     /* Computed Parameter: IOTypeSel_LoadIn_P1_l
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P2_Size_o[2];/* Computed Parameter: IOTypeSel_LoadIn_P2_Size_o
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P2_j;        /* Expression: FcnNos
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P3_Size_m[2];/* Computed Parameter: IOTypeSel_LoadIn_P3_Size_m
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P3_k;        /* Expression: width
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P4_Size_e[2];/* Computed Parameter: IOTypeSel_LoadIn_P4_Size_e
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T IOTypeSel_LoadIn_P4_o;        /* Expression: portType
                                        * Referenced by: '<S14>/IOTypeSel_LoadIn'
                                        */
  real_T Memory_X0_ko;                 /* Expression: 1
                                        * Referenced by: '<S14>/Memory'
                                        */
  real_T integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S38>/integrator'
                                        */
  real_T TransportDelay_Delay;         /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S38>/Transport Delay'
                                        */
  real_T K1_Value;                     /* Expression: 1./sps.Freq
                                        * Referenced by: '<S38>/K1'
                                        */
  real_T Memory_X0_a;                  /* Expression: sps.Vinit
                                        * Referenced by: '<S38>/Memory'
                                        */
  real_T integrator_IC_g;              /* Expression: 0
                                        * Referenced by: '<S39>/integrator'
                                        */
  real_T TransportDelay_Delay_b;       /* Expression: 1./sps.Freq
                                        * Referenced by: '<S39>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput_f;  /* Expression: 0
                                        * Referenced by: '<S39>/Transport Delay'
                                        */
  real_T K1_Value_k;                   /* Expression: 1./sps.Freq
                                        * Referenced by: '<S39>/K1'
                                        */
  real_T Memory_X0_c;                  /* Expression: sps.Vinit
                                        * Referenced by: '<S39>/Memory'
                                        */
  real_T OpWriteFile_P1_Size[2];       /* Computed Parameter: OpWriteFile_P1_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P1;               /* Expression: Acq_Group
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P2_Size[2];       /* Computed Parameter: OpWriteFile_P2_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P2[17];           /* Computed Parameter: OpWriteFile_P2
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P3_Size[2];       /* Computed Parameter: OpWriteFile_P3_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P3[5];            /* Computed Parameter: OpWriteFile_P3
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P4_Size[2];       /* Computed Parameter: OpWriteFile_P4_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P4;               /* Expression: Decimation
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P5_Size[2];       /* Computed Parameter: OpWriteFile_P5_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P5;               /* Expression: Nb_Samples
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P6_Size[2];       /* Computed Parameter: OpWriteFile_P6_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P6;               /* Expression: Buffer_size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P7_Size[2];       /* Computed Parameter: OpWriteFile_P7_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P7;               /* Expression: Sim_Mode
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P8_Size[2];       /* Computed Parameter: OpWriteFile_P8_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P8;               /* Expression: Static_File
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P9_Size[2];       /* Computed Parameter: OpWriteFile_P9_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P9;               /* Expression: file_size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P10_Size[2];      /* Computed Parameter: OpWriteFile_P10_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P10;              /* Expression: write_offline
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P11_Size[2];      /* Computed Parameter: OpWriteFile_P11_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P11[2];           /* Computed Parameter: OpWriteFile_P11
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P12_Size[2];      /* Computed Parameter: OpWriteFile_P12_Size
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T OpWriteFile_P12[2];           /* Computed Parameter: OpWriteFile_P12
                                        * Referenced by: '<S2>/OpWriteFile'
                                        */
  real_T Step_Time;                    /* Expression: 0.05
                                        * Referenced by: '<S12>/Step'
                                        */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S12>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 1
                                        * Referenced by: '<S12>/Step'
                                        */
  real_T RTEConversion_P1_Size[2];     /* Computed Parameter: RTEConversion_P1_Size
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P1;             /* Expression: nbMaxEvents
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P2_Size[2];     /* Computed Parameter: RTEConversion_P2_Size
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P2;             /* Expression: inputdatatype
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P3_Size[2];     /* Computed Parameter: RTEConversion_P3_Size
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P3;             /* Expression: outputdatatype
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P4_Size[2];     /* Computed Parameter: RTEConversion_P4_Size
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P4;             /* Expression: compensation
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P5_Size[2];     /* Computed Parameter: RTEConversion_P5_Size
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T RTEConversion_P5;             /* Expression: sampleTime
                                        * Referenced by: '<S12>/RTE Conversion'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0.1684
                                        * Referenced by: '<S62>/Unit Delay'
                                        */
  real_T Constant3_Value;              /* Expression: sps.Delay
                                        * Referenced by: '<S140>/Constant3'
                                        */
  real_T Constant1_Value_l;            /* Expression: sps.Period
                                        * Referenced by: '<S140>/Constant1'
                                        */
  real_T ib1_Gain;                     /* Expression: sps.Freq
                                        * Referenced by: '<S140>/1\ib1'
                                        */
  real_T LookupTable_XData[3];         /* Expression: [0 .5 1]
                                        * Referenced by: '<S140>/Lookup Table'
                                        */
  real_T LookupTable_YData[3];         /* Expression: [0 2 0]
                                        * Referenced by: '<S140>/Lookup Table'
                                        */
  real_T Constant2_Value_f;            /* Expression: 1
                                        * Referenced by: '<S140>/Constant2'
                                        */
  real_T Gain1_Gain_a;                 /* Expression: 0.5
                                        * Referenced by: '<S118>/Gain1'
                                        */
  real_T RTEConversion1_P1_Size[2];    /* Computed Parameter: RTEConversion1_P1_Size
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P1;            /* Expression: nbMaxEvents
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P2_Size[2];    /* Computed Parameter: RTEConversion1_P2_Size
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P2;            /* Expression: inputdatatype
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P3_Size[2];    /* Computed Parameter: RTEConversion1_P3_Size
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P3;            /* Expression: outputdatatype
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P4_Size[2];    /* Computed Parameter: RTEConversion1_P4_Size
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P4;            /* Expression: compensation
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P5_Size[2];    /* Computed Parameter: RTEConversion1_P5_Size
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTEConversion1_P5;            /* Expression: sampleTime
                                        * Referenced by: '<S12>/RTE Conversion1'
                                        */
  real_T RTELogicalOperator1_P1_Size[2];/* Computed Parameter: RTELogicalOperator1_P1_Size
                                         * Referenced by: '<S12>/RTE Logical Operator1'
                                         */
  real_T RTELogicalOperator1_P1;       /* Expression: LogicalOperator
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */
  real_T RTELogicalOperator1_P2_Size[2];/* Computed Parameter: RTELogicalOperator1_P2_Size
                                         * Referenced by: '<S12>/RTE Logical Operator1'
                                         */
  real_T RTELogicalOperator1_P2;       /* Expression: NbrInput
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */
  real_T RTELogicalOperator1_P3_Size[2];/* Computed Parameter: RTELogicalOperator1_P3_Size
                                         * Referenced by: '<S12>/RTE Logical Operator1'
                                         */
  real_T RTELogicalOperator1_P3;       /* Expression: NbrMaxEvents
                                        * Referenced by: '<S12>/RTE Logical Operator1'
                                        */
  real_T RTESPWM_P1_Size[2];           /* Computed Parameter: RTESPWM_P1_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P1;                   /* Expression: MaxEvents
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P2_Size[2];           /* Computed Parameter: RTESPWM_P2_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P2;                   /* Expression: SampleTime
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P3_Size[2];           /* Computed Parameter: RTESPWM_P3_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P3;                   /* Expression: MaxFreq
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P4_Size[2];           /* Computed Parameter: RTESPWM_P4_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P4;                   /* Expression: MinFreq
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P5_Size[2];           /* Computed Parameter: RTESPWM_P5_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P5;                   /* Expression: EnablingPort
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P6_Size[2];           /* Computed Parameter: RTESPWM_P6_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P6;                   /* Expression: NumberPhases
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P7_Size[2];           /* Computed Parameter: RTESPWM_P7_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P7;                   /* Expression: ComplementaryMode
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P8_Size[2];           /* Computed Parameter: RTESPWM_P8_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P8;                   /* Expression: RiseTimeDelay
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P9_Size[2];           /* Computed Parameter: RTESPWM_P9_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P9;                   /* Expression: CenterAlignmentMode
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P10_Size[2];          /* Computed Parameter: RTESPWM_P10_Size
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T RTESPWM_P10;                  /* Expression: SpaceVector
                                        * Referenced by: '<S12>/RTE SPWM'
                                        */
  real_T Switch_Threshold;             /* Expression: 1
                                        * Referenced by: '<S12>/Switch'
                                        */
  real_T RTEGround_P1_Size[2];         /* Computed Parameter: RTEGround_P1_Size
                                        * Referenced by: '<S12>/RTE Ground'
                                        */
  real_T RTEGround_P1;                 /* Expression: SampleTime
                                        * Referenced by: '<S12>/RTE Ground'
                                        */
  real_T RTE_Conversion_1_P1_Size[2];  /* Computed Parameter: RTE_Conversion_1_P1_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P1;          /* Expression: conversionType
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P2_Size[2];  /* Computed Parameter: RTE_Conversion_1_P2_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P2;          /* Expression: nbLine
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P3_Size[2];  /* Computed Parameter: RTE_Conversion_1_P3_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P3;          /* Expression: nbEvents
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P4_Size[2];  /* Computed Parameter: RTE_Conversion_1_P4_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P4;          /* Expression: timeUnit
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P5_Size[2];  /* Computed Parameter: RTE_Conversion_1_P5_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P5;          /* Expression: enTimeFactor
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P6_Size[2];  /* Computed Parameter: RTE_Conversion_1_P6_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P6;          /* Expression: sampleTime
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P7_Size[2];  /* Computed Parameter: RTE_Conversion_1_P7_Size
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T RTE_Conversion_1_P7;          /* Expression: initStates
                                        * Referenced by: '<S49>/RTE_Conversion_1'
                                        */
  real_T EventGen_eHS_1_P1_Size[2];    /* Computed Parameter: EventGen_eHS_1_P1_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P1[6];         /* Computed Parameter: EventGen_eHS_1_P1
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P2_Size[2];    /* Computed Parameter: EventGen_eHS_1_P2_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P2;            /* Expression: portNb
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P3_Size[2];    /* Computed Parameter: EventGen_eHS_1_P3_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P4_Size[2];    /* Computed Parameter: EventGen_eHS_1_P4_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P4;            /* Expression: size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P5_Size[2];    /* Computed Parameter: EventGen_eHS_1_P5_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P5;            /* Expression: nbChannels
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P6_Size[2];    /* Computed Parameter: EventGen_eHS_1_P6_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P6;            /* Expression: numwidth
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P7_Size[2];    /* Computed Parameter: EventGen_eHS_1_P7_Size
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T EventGen_eHS_1_P7;            /* Expression: timeUnit
                                        * Referenced by: '<S49>/EventGen_eHS_1'
                                        */
  real_T Constant1_Value_p;            /* Expression: 8.55
                                        * Referenced by: '<S12>/Constant1'
                                        */
  real_T Constant2_Value_c;            /* Expression: 37.4*14
                                        * Referenced by: '<S12>/Constant2'
                                        */
  real_T Constant3_Value_b;            /* Expression: 0.06
                                        * Referenced by: '<S12>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: 60*14
                                        * Referenced by: '<S12>/Constant4'
                                        */
  real_T UnitDelay_InitialCondition_a; /* Expression: sps.Finit
                                        * Referenced by: '<S87>/Unit Delay'
                                        */
  real_T avoiddivisionbyzero_UpperSat; /* Expression: 70
                                        * Referenced by: '<S76>/avoid division by zero'
                                        */
  real_T avoiddivisionbyzero_LowerSat; /* Expression: 40
                                        * Referenced by: '<S76>/avoid division by zero'
                                        */
  real_T Gain_Gain_e;                  /* Expression: 1/4
                                        * Referenced by: '<S76>/Gain'
                                        */
  real_T SFunction_P1_Size_i[2];       /* Computed Parameter: SFunction_P1_Size_i
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P1_ki;              /* Expression: MaxDelay
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P2_Size_i[2];       /* Computed Parameter: SFunction_P2_Size_i
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P2_k;               /* Expression: Ts
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P3_Size_i[2];       /* Computed Parameter: SFunction_P3_Size_i
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P3_n;               /* Expression: InitialValue
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P4_Size_n[2];       /* Computed Parameter: SFunction_P4_Size_n
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T SFunction_P4_a;               /* Expression: DFT
                                        * Referenced by: '<S113>/S-Function'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S87>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S87>/Discrete-Time Integrator'
                                        */
  real_T Constant4_Value_o;            /* Expression: 2*pi
                                        * Referenced by: '<S87>/Constant4'
                                        */
  real_T FirstcycleofsimulationId092Iq0_;/* Expression: 0
                                          * Referenced by: '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T FirstcycleofsimulationId092Iq_k;/* Expression: 1
                                          * Referenced by: '<S76>/First cycle of simulation Id=0.92, Iq=0'
                                          */
  real_T Duk_Gain;                     /* Expression: sps.D
                                        * Referenced by: '<S108>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition[2]; /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S108>/Delay_x1'
                                        */
  real_T C11_Gain;                     /* Expression: sps.C11
                                        * Referenced by: '<S111>/C11'
                                        */
  real_T Delay_x2_InitialCondition[2]; /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S108>/Delay_x2'
                                        */
  real_T C12_Gain;                     /* Expression: sps.C12
                                        * Referenced by: '<S111>/C12'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S62>/Unit Delay2'
                                        */
  real_T Iq_ref_Value;                 /* Expression: 0
                                        * Referenced by: '<S62>/Iq_ref'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S71>/Integrator'
                                        */
  real_T Integrator_IC;                /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S71>/Integrator'
                                        */
  real_T Gain1_Gain_b2;                /* Expression: 2
                                        * Referenced by: '<S72>/Gain1'
                                        */
  real_T Integ4_gainval_l;             /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S79>/Integ4'
                                        */
  real_T Integ4_IC_e;                  /* Expression: 0
                                        * Referenced by: '<S79>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_d;/* Expression: 1e6
                                          * Referenced by: '<S79>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_o;/* Expression: eps
                                          * Referenced by: '<S79>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_l[2];       /* Computed Parameter: SFunction_P1_Size_l
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P1_ke;              /* Expression: MaxDelay
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P2_Size_c[2];       /* Computed Parameter: SFunction_P2_Size_c
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P2_db;              /* Expression: Ts
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P3_Size_b[2];       /* Computed Parameter: SFunction_P3_Size_b
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P3_m;               /* Expression: InitialValue
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P4_Size_f[2];       /* Computed Parameter: SFunction_P4_Size_f
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T SFunction_P4_aj;              /* Expression: DFT
                                        * Referenced by: '<S81>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_f; /* Expression: 0
                                        * Referenced by: '<S80>/Unit Delay'
                                        */
  real_T Constant_Value_p;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S79>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_j;/* Expression: sps.Vinit
                                        * Referenced by: '<S79>/Unit Delay1'
                                        */
  real_T Gain3_Gain_k;                 /* Expression: 2
                                        * Referenced by: '<S72>/Gain3'
                                        */
  real_T Integ4_gainval_g;             /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S82>/Integ4'
                                        */
  real_T Integ4_IC_l;                  /* Expression: 0
                                        * Referenced by: '<S82>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_e;/* Expression: 1e6
                                          * Referenced by: '<S82>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerS_h5;/* Expression: eps
                                          * Referenced by: '<S82>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_hi[2];      /* Computed Parameter: SFunction_P1_Size_hi
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P1_ks;              /* Expression: MaxDelay
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P2_Size_j[2];       /* Computed Parameter: SFunction_P2_Size_j
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P2_e;               /* Expression: Ts
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P3_Size_c[2];       /* Computed Parameter: SFunction_P3_Size_c
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P3_c;               /* Expression: InitialValue
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P4_Size_m[2];       /* Computed Parameter: SFunction_P4_Size_m
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T SFunction_P4_g;               /* Expression: DFT
                                        * Referenced by: '<S84>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_d; /* Expression: 0
                                        * Referenced by: '<S83>/Unit Delay'
                                        */
  real_T Constant_Value_d;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S82>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_e;/* Expression: sps.Vinit
                                        * Referenced by: '<S82>/Unit Delay1'
                                        */
  real_T RadDeg_Gain_n;                /* Expression: 180/pi
                                        * Referenced by: '<S72>/Rad->Deg.'
                                        */
  real_T torad_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S67>/to-rad'
                                        */
  real_T Rff_Gain;                     /* Expression: RLff(1)
                                        * Referenced by: '<S64>/Rff '
                                        */
  real_T Lff_Gain;                     /* Expression: RLff(2)
                                        * Referenced by: '<S64>/Lff  '
                                        */
  real_T Rff_Gain_o;                   /* Expression: RLff(1)
                                        * Referenced by: '<S64>/Rff'
                                        */
  real_T Lff_Gain_b;                   /* Expression: RLff(2)
                                        * Referenced by: '<S64>/Lff'
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: 1.5
                                        * Referenced by: '<S64>/Saturation'
                                        */
  real_T Saturation_LowerSat_e;        /* Expression: -1.5
                                        * Referenced by: '<S64>/Saturation'
                                        */
  real_T Iph_1_Value;                  /* Expression: Limits_MPPT(1)
                                        * Referenced by: '<S65>/Iph_1'
                                        */
  real_T Iph_2_Value;                  /* Expression: Limits_MPPT(2)
                                        * Referenced by: '<S65>/Iph_2'
                                        */
  real_T MPPT_On_Value;                /* Expression: 1
                                        * Referenced by: '<S12>/MPPT_On'
                                        */
  real_T Integ4_gainval_o;             /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S85>/Integ4'
                                        */
  real_T Integ4_IC_o;                  /* Expression: 0
                                        * Referenced by: '<S85>/Integ4'
                                        */
  real_T K1_Value_f;                   /* Expression: sps.Delay
                                        * Referenced by: '<S85>/K1'
                                        */
  real_T SFunction_P1_Size_o[2];       /* Computed Parameter: SFunction_P1_Size_o
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P1_cs;              /* Expression: MaxDelay
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P2_Size_id[2];      /* Computed Parameter: SFunction_P2_Size_id
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P2_a;               /* Expression: Ts
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P3_Size_l[2];       /* Computed Parameter: SFunction_P3_Size_l
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P3_f;               /* Expression: InitialValue
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P4_Size_i[2];       /* Computed Parameter: SFunction_P4_Size_i
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T SFunction_P4_dl;              /* Expression: DFT
                                        * Referenced by: '<S86>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_m; /* Expression: 0
                                        * Referenced by: '<S85>/Unit Delay'
                                        */
  real_T K2_Value;                     /* Expression: sps.Freq
                                        * Referenced by: '<S85>/K2'
                                        */
  real_T UnitDelay1_InitialCondition_m;/* Expression: sps.Vinit
                                        * Referenced by: '<S85>/Unit Delay1'
                                        */
  real_T Integ4_gainval_b;             /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S101>/Integ4'
                                        */
  real_T Integ4_IC_m;                  /* Expression: 0
                                        * Referenced by: '<S101>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperS_ix;/* Expression: 1e6
                                          * Referenced by: '<S101>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_e;/* Expression: eps
                                          * Referenced by: '<S101>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_lc[2];      /* Computed Parameter: SFunction_P1_Size_lc
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P1_g;               /* Expression: MaxDelay
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P2_Size_g[2];       /* Computed Parameter: SFunction_P2_Size_g
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P2_ac;              /* Expression: Ts
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P3_Size_f[2];       /* Computed Parameter: SFunction_P3_Size_f
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P3_nw;              /* Expression: InitialValue
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P4_Size_a[2];       /* Computed Parameter: SFunction_P4_Size_a
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T SFunction_P4_d1;              /* Expression: DFT
                                        * Referenced by: '<S103>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_e; /* Expression: 0
                                        * Referenced by: '<S102>/Unit Delay'
                                        */
  real_T Constant_Value_c;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S101>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_m1;/* Expression: sps.Vinit
                                         * Referenced by: '<S101>/Unit Delay1'
                                         */
  real_T DiscreteDerivative_DenCoef[2];/* Expression: [ TcD  Ts-TcD ]
                                        * Referenced by: '<S89>/Discrete Derivative '
                                        */
  real_T DiscreteDerivative_InitialState;/* Expression: 0
                                          * Referenced by: '<S89>/Discrete Derivative '
                                          */
  real_T DiscreteTimeIntegrator_gainva_a;/* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                          * Referenced by: '<S89>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperSat;/* Expression: Par_Limits(1)
                                          * Referenced by: '<S89>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerSat;/* Expression: Par_Limits(2)
                                          * Referenced by: '<S89>/Discrete-Time Integrator'
                                          */
  real_T Saturation1_UpperSat_c;       /* Expression: Par_Limits(1)
                                        * Referenced by: '<S89>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_j;       /* Expression: Par_Limits(2)
                                        * Referenced by: '<S89>/Saturation1'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1/2/pi
                                        * Referenced by: '<S87>/Gain10'
                                        */
  real_T RateLimiter_RisingLim;        /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  real_T RateLimiter_IC;               /* Expression: sps.Finit
                                        * Referenced by: '<S87>/Rate Limiter'
                                        */
  real_T Delay_x1_InitialCondition_b;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S104>/Delay_x1'
                                        */
  real_T A11_Gain;                     /* Expression: sps.A11
                                        * Referenced by: '<S105>/A11'
                                        */
  real_T Delay_x2_InitialCondition_h;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S104>/Delay_x2'
                                        */
  real_T A12_Gain;                     /* Expression: sps.A12
                                        * Referenced by: '<S105>/A12'
                                        */
  real_T A21_Gain;                     /* Expression: sps.A21
                                        * Referenced by: '<S105>/A21'
                                        */
  real_T A22_Gain;                     /* Expression: sps.A22
                                        * Referenced by: '<S105>/A22'
                                        */
  real_T B11_Gain;                     /* Expression: sps.B11
                                        * Referenced by: '<S106>/B11'
                                        */
  real_T B21_Gain;                     /* Expression: sps.B21
                                        * Referenced by: '<S106>/B21'
                                        */
  real_T Duk_Gain_b;                   /* Expression: sps.D
                                        * Referenced by: '<S104>/D*u(k)'
                                        */
  real_T C11_Gain_b;                   /* Expression: sps.C11
                                        * Referenced by: '<S107>/C11'
                                        */
  real_T C12_Gain_f;                   /* Expression: sps.C12
                                        * Referenced by: '<S107>/C12'
                                        */
  real_T A11_Gain_k;                   /* Expression: sps.A11
                                        * Referenced by: '<S109>/A11'
                                        */
  real_T A12_Gain_k;                   /* Expression: sps.A12
                                        * Referenced by: '<S109>/A12'
                                        */
  real_T A21_Gain_d;                   /* Expression: sps.A21
                                        * Referenced by: '<S109>/A21'
                                        */
  real_T A22_Gain_k;                   /* Expression: sps.A22
                                        * Referenced by: '<S109>/A22'
                                        */
  real_T B11_Gain_f;                   /* Expression: sps.B11
                                        * Referenced by: '<S110>/B11'
                                        */
  real_T B21_Gain_j;                   /* Expression: sps.B21
                                        * Referenced by: '<S110>/B21'
                                        */
  real_T Constant1_Value_a;            /* Expression: 0
                                        * Referenced by: '<S76>/Constant1'
                                        */
  real_T Constant2_Value_i;            /* Expression: 0
                                        * Referenced by: '<S69>/Constant2'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S62>/Unit Delay3'
                                        */
  real_T Gain1_Gain_ne;                /* Expression: 1
                                        * Referenced by: '<S69>/Gain1'
                                        */
  real_T UnitDelay1_InitialCondition_eb;/* Expression: 0
                                         * Referenced by: '<S62>/Unit Delay1'
                                         */
  real_T Integrator_gainval_e;         /* Computed Parameter: Integrator_gainval_e
                                        * Referenced by: '<S142>/Integrator'
                                        */
  real_T Integrator_IC_p;              /* Expression: InitialConditionForIntegrator
                                        * Referenced by: '<S142>/Integrator'
                                        */
  real_T RTEPeriodMeter_P1_Size[2];    /* Computed Parameter: RTEPeriodMeter_P1_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P1;            /* Expression: Period
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P2_Size[2];    /* Computed Parameter: RTEPeriodMeter_P2_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P2;            /* Expression: Frequency
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P3_Size[2];    /* Computed Parameter: RTEPeriodMeter_P3_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P3;            /* Expression: EdgeType
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P4_Size[2];    /* Computed Parameter: RTEPeriodMeter_P4_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P4;            /* Expression: TOn
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P5_Size[2];    /* Computed Parameter: RTEPeriodMeter_P5_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P5;            /* Expression: TOff
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P6_Size[2];    /* Computed Parameter: RTEPeriodMeter_P6_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P6;            /* Expression: DutyCycle
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P7_Size[2];    /* Computed Parameter: RTEPeriodMeter_P7_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P7;            /* Expression: DuringType
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P8_Size[2];    /* Computed Parameter: RTEPeriodMeter_P8_Size
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T RTEPeriodMeter_P8;            /* Expression: MinimumFrequency
                                        * Referenced by: '<S12>/RTE Period Meter'
                                        */
  real_T Saturation_UpperSat_b1;       /* Expression: 10
                                        * Referenced by: '<S12>/Saturation'
                                        */
  real_T Saturation_LowerSat_c;        /* Expression: 0
                                        * Referenced by: '<S12>/Saturation'
                                        */
  real_T Constant1_Value_h;            /* Expression: 0.5
                                        * Referenced by: '<S2>/Constant1'
                                        */
  real_T OpTrigger_P1_Size[2];         /* Computed Parameter: OpTrigger_P1_Size
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpTrigger_P1;                 /* Expression: Acq_Group
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpTrigger_P2_Size[2];         /* Computed Parameter: OpTrigger_P2_Size
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpTrigger_P2;                 /* Expression: Trig_Type
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpTrigger_P3_Size[2];         /* Computed Parameter: OpTrigger_P3_Size
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpTrigger_P3;                 /* Expression: Trig_Offset
                                        * Referenced by: '<S2>/OpTrigger'
                                        */
  real_T OpCtrl_P1_Size[2];            /* Computed Parameter: OpCtrl_P1_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P1[6];                 /* Computed Parameter: OpCtrl_P1
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P2_Size[2];            /* Computed Parameter: OpCtrl_P2_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P2;                    /* Expression: boardid
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P3_Size[2];            /* Computed Parameter: OpCtrl_P3_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P3;                    /* Expression: mode
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P4_Size[2];            /* Computed Parameter: OpCtrl_P4_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P4;                    /* Expression: externalClock
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P5_Size[2];            /* Computed Parameter: OpCtrl_P5_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P5;                    /* Expression: decimRtsi
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P6_Size[2];            /* Computed Parameter: OpCtrl_P6_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P6;                    /* Expression: 1
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P7_Size[2];            /* Computed Parameter: OpCtrl_P7_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P7;                    /* Expression: SampleTime
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P8_Size[2];            /* Computed Parameter: OpCtrl_P8_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P8;                    /* Expression: calibIO
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P9_Size[2];            /* Computed Parameter: OpCtrl_P9_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P9;                    /* Expression: numconfig
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P10_Size[2];           /* Computed Parameter: OpCtrl_P10_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P10;                   /* Expression: loadinport
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P11_Size[2];           /* Computed Parameter: OpCtrl_P11_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P11;                   /* Expression: BoardType
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P12_Size[2];           /* Computed Parameter: OpCtrl_P12_Size
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  real_T OpCtrl_P12;                   /* Expression: sync_type
                                        * Referenced by: '<S2>/OpCtrl'
                                        */
  uint32_T default_Value[17];          /* Computed Parameter: default_Value
                                        * Referenced by: '<S40>/default'
                                        */
  uint32_T addr4_Value;                /* Computed Parameter: addr4_Value
                                        * Referenced by: '<S40>/addr4'
                                        */
  uint32_T addr3_Value;                /* Computed Parameter: addr3_Value
                                        * Referenced by: '<S40>/addr3'
                                        */
  uint32_T addr2_Value;                /* Computed Parameter: addr2_Value
                                        * Referenced by: '<S40>/addr2'
                                        */
  uint32_T addr1_Value;                /* Computed Parameter: addr1_Value
                                        * Referenced by: '<S40>/addr1'
                                        */
  uint32_T addr_Value;                 /* Computed Parameter: addr_Value
                                        * Referenced by: '<S40>/addr'
                                        */
  uint32_T sat_scn_UpperSat;           /* Computed Parameter: sat_scn_UpperSat
                                        * Referenced by: '<S10>/sat_scn'
                                        */
  uint32_T sat_scn_LowerSat;           /* Computed Parameter: sat_scn_LowerSat
                                        * Referenced by: '<S10>/sat_scn'
                                        */
  uint32_T blockID_Value;              /* Computed Parameter: blockID_Value
                                        * Referenced by: '<S40>/blockID'
                                        */
  uint32_T IOTypeSel_Value;            /* Computed Parameter: IOTypeSel_Value
                                        * Referenced by: '<S13>/IOTypeSel'
                                        */
  uint32_T IOTypeSel_Value_n;          /* Computed Parameter: IOTypeSel_Value_n
                                        * Referenced by: '<S14>/IOTypeSel'
                                        */
  uint32_T Saturation_UpperSat_k;      /* Computed Parameter: Saturation_UpperSat_k
                                        * Referenced by: '<S10>/Saturation'
                                        */
  uint32_T Saturation_LowerSat_f;      /* Computed Parameter: Saturation_LowerSat_f
                                        * Referenced by: '<S10>/Saturation'
                                        */
  uint32_T load_config1_Value[10];     /* Computed Parameter: load_config1_Value
                                        * Referenced by: '<S10>/load_config1'
                                        */
  uint32_T shift_2bits_Gain;           /* Computed Parameter: shift_2bits_Gain
                                        * Referenced by: '<S10>/shift_2bits'
                                        */
  boolean_T Constant_Value_gg;         /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S45>/Constant'
                                        */
  boolean_T Constant_Value_n;          /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S43>/Constant'
                                        */
  boolean_T Constant_Value_kj;         /* Expression: ~strcmp(FormatType,'Unsigned')
                                        * Referenced by: '<S44>/Constant'
                                        */
  P_TrueRMS_boost_and_two_level_T TrueRMS_i;/* '<S7>/TrueRMS ' */
  P_RMS_boost_and_two_level__1__T RMS_h;/* '<S7>/RMS ' */
  P_TrueRMS_boost_and_two_level_T TrueRMS;/* '<S6>/TrueRMS ' */
  P_RMS_boost_and_two_level__1__T RMS; /* '<S6>/RMS ' */
};

/* Real-time Model Data Structure */
struct tag_RTM_boost_and_two_level__1_sm_ehs_T {
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
    SimStruct childSFunctions[42];
    SimStruct *childSFunctionPtrs[42];
    struct _ssBlkInfo2 blkInfo2[42];
    struct _ssSFcnModelMethods2 methods2[42];
    struct _ssSFcnModelMethods3 methods3[42];
    struct _ssStatesInfo2 statesInfo2[42];
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
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn2;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[2];
      mxArray *params[2];
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn4;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[10];
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
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[3];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn6;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[8];
      mxArray *params[8];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn7;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[2];
      uint_T attribs[9];
      mxArray *params[9];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn8;

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
    } Sfcn9;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[8];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn10;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[3];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn11;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[6];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn12;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[2];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn13;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[16];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn14;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[16];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn15;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn16;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[2];
      uint_T attribs[9];
      mxArray *params[9];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn17;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[7];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn18;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[16];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn19;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[8];
      mxArray *params[8];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn20;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[3];
      uint_T attribs[6];
      mxArray *params[6];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn21;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn22;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[3];
      uint_T attribs[6];
      mxArray *params[6];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn23;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn24;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[7];
      struct _ssPortOutputs outputPortInfo[1];
    } Sfcn25;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[6];
      struct _ssPortOutputs outputPortInfo[2];
      uint_T attribs[12];
      mxArray *params[12];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn26;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[5];
      mxArray *params[5];
      struct _ssDWorkRecord dWork[2];
      struct _ssDWorkAuxRecord dWorkAux[2];
    } Sfcn27;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[5];
      mxArray *params[5];
      struct _ssDWorkRecord dWork[2];
      struct _ssDWorkAuxRecord dWorkAux[2];
    } Sfcn28;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[4];
      real_T const *UPtrs1[4];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn29;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      struct _ssPortOutputs outputPortInfo[2];
      uint_T attribs[10];
      mxArray *params[10];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn30;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn31;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[8];
      struct _ssPortOutputs outputPortInfo[16];
      uint_T attribs[7];
      mxArray *params[7];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn32;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[16];
      real_T const *UPtrs0[4];
      real_T const *UPtrs1[4];
      real_T const *UPtrs2[4];
      real_T const *UPtrs3[4];
      real_T const *UPtrs4[4];
      real_T const *UPtrs5[4];
      real_T const *UPtrs6[4];
      real_T const *UPtrs7[4];
      real_T const *UPtrs8[4];
      real_T const *UPtrs9[4];
      real_T const *UPtrs10[4];
      real_T const *UPtrs11[4];
      real_T const *UPtrs12[4];
      real_T const *UPtrs13[4];
      real_T const *UPtrs14[4];
      real_T const *UPtrs15[4];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[7];
      mxArray *params[7];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn33;

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
    } Sfcn34;

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
    } Sfcn35;

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
    } Sfcn36;

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
    } Sfcn37;

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
    } Sfcn38;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[8];
      mxArray *params[8];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn39;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[3];
      mxArray *params[3];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn40;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[2];
      uint_T attribs[12];
      mxArray *params[12];
      struct _ssDWorkRecord dWork[2];
      struct _ssDWorkAuxRecord dWorkAux[2];
    } Sfcn41;
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
    real_T odeY[8];
    real_T odeF[3][8];
    ODE3_IntgData intgData;
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
extern P_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_P;

/* Block signals (auto storage) */
extern B_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_B;

/* Continuous states (auto storage) */
extern X_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_X;

/* Block states (auto storage) */
extern DW_boost_and_two_level__1_sm_ehs_T boost_and_two_level__1_sm_ehs_DW;

/* External data declarations for dependent source files */
extern const real_T boost_and_two_level__1_sm_ehs_RGND;/* real_T ground */

/* Model entry point functions */
extern void boost_and_two_level__1_sm_ehs_initialize(void);
extern void boost_and_two_level__1_sm_ehs_output(void);
extern void boost_and_two_level__1_sm_ehs_update(void);
extern void boost_and_two_level__1_sm_ehs_terminate(void);

/*====================*
 * External functions *
 *====================*/
extern boost_and_two_level__1_sm_ehs_rtModel *boost_and_two_level__1_sm_ehs(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_boost_and_two_level__1_sm_ehs_T *const
  boost_and_two_level__1_sm_ehs_M;

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
 * '<Root>' : 'boost_and_two_level__1_sm_ehs'
 * '<S1>'   : 'boost_and_two_level__1_sm_ehs/OpCCode_do_not_touch'
 * '<S2>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS'
 * '<S3>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/Compare To Constant'
 * '<S4>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/From Digital Inputs (Pulse-Width Analyzers)'
 * '<S5>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/OpComm'
 * '<S6>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS'
 * '<S7>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1'
 * '<S8>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/Subsystem'
 * '<S9>'   : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs'
 * '<S10>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk'
 * '<S11>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem'
 * '<S12>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates'
 * '<S13>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/From Digital Inputs (Pulse-Width Analyzers)/Selectable DI'
 * '<S14>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/From Digital Inputs (Pulse-Width Analyzers)/Selectable DI1'
 * '<S15>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/OpComm/Receive'
 * '<S16>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/OpComm/busStruct'
 * '<S17>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/OpComm/busStruct/Sub1'
 * '<S18>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS '
 * '<S19>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/TrueRMS '
 * '<S20>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS /Fourier1'
 * '<S21>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS /Fourier1/Mean'
 * '<S22>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS /Fourier1/Mean value1'
 * '<S23>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS /Fourier1/Mean/Model'
 * '<S24>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/RMS /Fourier1/Mean value1/Model'
 * '<S25>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/TrueRMS /Mean value'
 * '<S26>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS/TrueRMS /Mean value/Model'
 * '<S27>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS '
 * '<S28>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/TrueRMS '
 * '<S29>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS /Fourier1'
 * '<S30>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS /Fourier1/Mean'
 * '<S31>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS /Fourier1/Mean value1'
 * '<S32>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS /Fourier1/Mean/Model'
 * '<S33>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/RMS /Fourier1/Mean value1/Model'
 * '<S34>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/TrueRMS /Mean value'
 * '<S35>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/RMS1/TrueRMS /Mean value/Model'
 * '<S36>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/Subsystem/Mean'
 * '<S37>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/Subsystem/Mean1'
 * '<S38>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/Subsystem/Mean/Model'
 * '<S39>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/Subsystem/Mean1/Model'
 * '<S40>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments'
 * '<S41>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Subsystem'
 * '<S42>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Counter'
 * '<S43>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Rescale to  Fixed-Pt  format'
 * '<S44>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Rescale to  Fixed-Pt  format1'
 * '<S45>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/Rescale to  Fixed-Pt  format2'
 * '<S46>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/change_detector'
 * '<S47>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/To Analog Outputs/Aout_Mapping&Adjustments/change_detector/Detect Change'
 * '<S48>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/Detect Change'
 * '<S49>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/c_solver_rte1'
 * '<S50>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/change_detector'
 * '<S51>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/eHSx64 Gen3 CommBlk/change_detector/Detect Change'
 * '<S52>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem1'
 * '<S53>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem2'
 * '<S54>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem3'
 * '<S55>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem4'
 * '<S56>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem5'
 * '<S57>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem1/Send1'
 * '<S58>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem2/Send2'
 * '<S59>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem3/Send3'
 * '<S60>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem4/Send4'
 * '<S61>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/rtlab_send_subsystem/Subsystem5/Send5'
 * '<S62>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control'
 * '<S63>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/MATLAB Function'
 * '<S64>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/Current Regulator'
 * '<S65>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/MPPT  Parameters'
 * '<S66>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/MPPT Controller using Perturbe  & Observe technique  '
 * '<S67>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements'
 * '<S68>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator'
 * '<S69>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/U_ref Generation '
 * '<S70>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/VDC Regulator'
 * '<S71>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/Current Regulator/PI'
 * '<S72>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)'
 * '<S73>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Mean'
 * '<S74>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL'
 * '<S75>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Second-Order Filter'
 * '<S76>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform'
 * '<S77>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S78>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S79>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S80>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S81>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S82>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S83>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S84>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S85>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Mean/Model'
 * '<S86>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Mean/Model/Discrete Variable Time Delay'
 * '<S87>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model'
 * '<S88>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control'
 * '<S89>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Discrete'
 * '<S90>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)'
 * '<S91>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter'
 * '<S92>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)'
 * '<S93>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S94>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S95>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S96>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S97>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S98>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S99>'  : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S100>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Automatic Gain Control/Fundamental (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S101>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model'
 * '<S102>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Correction subsystem'
 * '<S103>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Mean (Variable Frequency)/Model/Discrete Variable Time Delay'
 * '<S104>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model'
 * '<S105>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/A*k(k-1)'
 * '<S106>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S107>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/PLL/Model/Second-Order Filter/Model/C*x(k)'
 * '<S108>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Second-Order Filter/Model'
 * '<S109>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Second-Order Filter/Model/A*k(k-1)'
 * '<S110>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S111>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Second-Order Filter/Model/C*x(k)'
 * '<S112>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0'
 * '<S113>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Discrete Variable Time Delay'
 * '<S114>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S115>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S116>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S117>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PLL & Measurements/Single-Phase dq Transform/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S118>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Cr_MinMax'
 * '<S119>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Modulator type'
 * '<S120>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Reference signal'
 * '<S121>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling'
 * '<S122>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Modulator type/Full Bridge Bipolar'
 * '<S123>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Modulator type/Full Bridge Unipolar'
 * '<S124>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Modulator type/One Three Phase Bridge'
 * '<S125>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Reference signal/External'
 * '<S126>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Reference signal/Internal'
 * '<S127>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical'
 * '<S128>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Natural'
 * '<S129>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical'
 * '<S130>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical'
 * '<S131>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Natural'
 * '<S132>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical'
 * '<S133>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Asymmetrical/Sample & Hold'
 * '<S134>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Natural/Sync_NaturalSampling'
 * '<S135>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling'
 * '<S136>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Sync Symmetrical/Sync_SymmetricalSampling/Sample & Hold'
 * '<S137>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Asymmetrical/Unsync_AsymmetricalSampling'
 * '<S138>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling'
 * '<S139>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator'
 * '<S140>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator/Model'
 * '<S141>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/PWM_Generator/Sampling/Unsync Symmetrical/Unsync_SymmetricalSampling'
 * '<S142>' : 'boost_and_two_level__1_sm_ehs/SM_eHS/source_and_gates/Inverter Control/VDC Regulator/PI'
 */
#endif                                 /* RTW_HEADER_boost_and_two_level__1_sm_ehs_h_ */
