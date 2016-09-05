/**
 * This function is generated by RT-LAB during model compilation (at 'generation' step).
 * This function copies data from src to dst, ignoring pointers.
 * Note that only PWork are supposed to be pointers within a DWork structure
 * and that sub-structures are copied in one operation since they do not contain pointers.
 */
int OpalSnapshot_Copy_DWork(void * src, void * dst) {
   D_Work * pSrc = (D_Work*)src;
   D_Work * pDst = (D_Work*)dst;
   int size = 0, eltSize = 0;

   eltSize = sizeof(pSrc->Integ4_DSTATE);
   memcpy(&pDst->Integ4_DSTATE, &pSrc->Integ4_DSTATE, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->UnitDelay_DSTATE);
   memcpy(&pDst->UnitDelay_DSTATE, &pSrc->UnitDelay_DSTATE, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->UnitDelay1_DSTATE);
   memcpy(&pDst->UnitDelay1_DSTATE, &pSrc->UnitDelay1_DSTATE, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Integ4_DSTATE_e);
   memcpy(&pDst->Integ4_DSTATE_e, &pSrc->Integ4_DSTATE_e, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->UnitDelay_DSTATE_l);
   memcpy(&pDst->UnitDelay_DSTATE_l, &pSrc->UnitDelay_DSTATE_l, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->UnitDelay1_DSTATE_p);
   memcpy(&pDst->UnitDelay1_DSTATE_p, &pSrc->UnitDelay1_DSTATE_p, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_PreviousInput);
   memcpy(&pDst->SFunction_PreviousInput, &pSrc->SFunction_PreviousInput, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Memory_PreviousInput);
   memcpy(&pDst->Memory_PreviousInput, &pSrc->Memory_PreviousInput, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_RWORK);
   memcpy(&pDst->SFunction_RWORK, &pSrc->SFunction_RWORK, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_RWORK_f);
   memcpy(&pDst->SFunction_RWORK_f, &pSrc->SFunction_RWORK_f, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_IWORK);
   memcpy(&pDst->SFunction_IWORK, &pSrc->SFunction_IWORK, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_IWORK_f);
   memcpy(&pDst->SFunction_IWORK_f, &pSrc->SFunction_IWORK_f, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->SFunction_IWORK_n);
   memcpy(&pDst->SFunction_IWORK_n, &pSrc->SFunction_IWORK_n, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Integ4_SYSTEM_ENABLE);
   memcpy(&pDst->Integ4_SYSTEM_ENABLE, &pSrc->Integ4_SYSTEM_ENABLE, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Integ4_SYSTEM_ENABLE_f);
   memcpy(&pDst->Integ4_SYSTEM_ENABLE_f, &pSrc->Integ4_SYSTEM_ENABLE_f, eltSize);
   size += eltSize;

   return size;
}

/**
 * This function is generated by RT-LAB during model compilation (at 'generation' step).
 * This function copies data from a raw buffer (src) 
 * to an RT_MODEL structure(dst), ignoring pointers.
 * Note that sub-structures are copied in one operation since they cannot contain pointers.
 */
int OpalSnapshot_Copy_TimingData(void * src, void * dst) {
   RT_MODEL	tmpBuffer;
   RT_MODEL * pSrc = &tmpBuffer;
   RT_MODEL * pDst = (RT_MODEL *)dst;
   int TimingSize = sizeof(tmpBuffer.Timing);
   int size = 0, eltSize = 0;

   memcpy(&tmpBuffer.Timing, src, TimingSize);

   eltSize = sizeof(pSrc->Timing.clockTick0);
   memcpy(&pDst->Timing.clockTick0, &pSrc->Timing.clockTick0, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.clockTickH0);
   memcpy(&pDst->Timing.clockTickH0, &pSrc->Timing.clockTickH0, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.tFinal);
   memcpy(&pDst->Timing.tFinal, &pSrc->Timing.tFinal, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.offsetTimesArray);
   memcpy(&pDst->Timing.offsetTimesArray, &pSrc->Timing.offsetTimesArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.perTaskSampleHitsArray);
   memcpy(&pDst->Timing.perTaskSampleHitsArray, &pSrc->Timing.perTaskSampleHitsArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.sampleHitArray);
   memcpy(&pDst->Timing.sampleHitArray, &pSrc->Timing.sampleHitArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.sampleTimesArray);
   memcpy(&pDst->Timing.sampleTimesArray, &pSrc->Timing.sampleTimesArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.sampleTimeTaskIDArray);
   memcpy(&pDst->Timing.sampleTimeTaskIDArray, &pSrc->Timing.sampleTimeTaskIDArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.simTimeStep);
   memcpy(&pDst->Timing.simTimeStep, &pSrc->Timing.simTimeStep, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.tStart);
   memcpy(&pDst->Timing.tStart, &pSrc->Timing.tStart, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.stepSize);
   memcpy(&pDst->Timing.stepSize, &pSrc->Timing.stepSize, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.stepSize0);
   memcpy(&pDst->Timing.stepSize0, &pSrc->Timing.stepSize0, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.stopRequestedFlag);
   memcpy(&pDst->Timing.stopRequestedFlag, &pSrc->Timing.stopRequestedFlag, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.tArray);
   memcpy(&pDst->Timing.tArray, &pSrc->Timing.tArray, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->Timing.timeOfLastOutput);
   memcpy(&pDst->Timing.timeOfLastOutput, &pSrc->Timing.timeOfLastOutput, eltSize);
   size += eltSize;

   return size;
}

/**
 * This function is generated by RT-LAB during model compilation (at 'generation' step).
 * This function copies data from a raw buffer (src) 
 * to an RT_MODEL structure(dst), ignoring pointers.
 * Note that ingData must not be copied since it contains pointers.
 */
int OpalSnapshot_Copy_ModelData(void * src, void * dst) {
   RT_MODEL	tmpBuffer;
   RT_MODEL * pSrc = &tmpBuffer;
   RT_MODEL * pDst = (RT_MODEL *)dst;
   int ModelDataSize = sizeof(tmpBuffer.ModelData);
   int size = 0, eltSize = 0;

   memcpy(&tmpBuffer.ModelData, src, ModelDataSize);

   eltSize = sizeof(pSrc->ModelData.blkStateChange);
   memcpy(&pDst->ModelData.blkStateChange, &pSrc->ModelData.blkStateChange, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->ModelData.derivCacheNeedsReset);
   memcpy(&pDst->ModelData.derivCacheNeedsReset, &pSrc->ModelData.derivCacheNeedsReset, eltSize);
   size += eltSize;

   eltSize = sizeof(pSrc->ModelData.zCCacheNeedsReset);
   memcpy(&pDst->ModelData.zCCacheNeedsReset, &pSrc->ModelData.zCCacheNeedsReset, eltSize);
   size += eltSize;

   return size;
}

