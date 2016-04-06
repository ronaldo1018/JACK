#include    "jack.h"
#include    <time.h>
#include    <malloc.h>
#include    <CL/cl.h>
#include    "cl_trace.h"
#include    "cl_util.h"
#include    "cl_set.h"

void JackCL(int NMC, int bat_pltsel, int bat_devsel, int bat_wgropitems, int bat_batchitems,
            REAL EMEAN, REAL ESTDV, REAL SIGMA, REAL XYANGLE, REAL AZIMUTH, collector *C,
            int *g_mem_size, int *c_mem_size, int *l_mem_size, 
            int *c_ProtonWater_Energy_size, int *g_traceA_size, int *g_arena_size)
{
//==============//  
// OpenCL begin //
//==============//
//========================//  
// Basic Platform Setting //
//========================//
  cl_int ret;
  size_t maxCmptUnits;
  cl_event rdEvent_trace, wrRanEvent, MkernelEvent;
  cl_ulong eventStart, eventEnd, totalTime=0;
  cl_ulong initime=0, kerneltime=0, rdbacktime=0;
  time_t t = time(0);
  srand48(t);
  
  
  if (bat_pltsel == 99)
  { 
      QueryHWinfo(&maxCmptUnits);  
      return;
  }
  else
  {
      GetHW();  
  }
// Input Working elements and relative data //

  int plt_sel, dev_sel;     
  size_t ItermsPerWGroup;
  size_t ItermsPerBatch;
  size_t Batches;

//  char *intrace = argv[1];
  //for batch File
//  if (argc >= 2 && intrace[0] == 'o')
//  {
      plt_sel = bat_pltsel;
      dev_sel = bat_devsel;
      ItermsPerWGroup = bat_wgropitems;
      ItermsPerBatch = bat_batchitems;
//  }
//  else
//  {
//      cl_OperateSelect(&plt_sel, &dev_sel, NMC, EMEAN, ESTDV, 
//                   &ItermsPerWGroup, &ItermsPerBatch);
//  }
  
  Batches = NMC/ItermsPerBatch;
  
//  printf("Starting Jack\n");
// Check and Create OpenCL Context & CmdQ //

  cl_context context;
  cl_command_queue command_queue;
  cl_command_queue cmdQ_Read;
  
  SetCont(plt_sel, dev_sel, &context);
  
  command_queue = clCreateCommandQueue(context, 
                  Devices[dev_sel], CL_QUEUE_PROFILING_ENABLE, &ret);
  cl_errChk(ret,"Error >> clCreateCommandQueue");

  cmdQ_Read = clCreateCommandQueue(context, 
            Devices[dev_sel], CL_QUEUE_PROFILING_ENABLE, &ret);
  cl_errChk(ret,"Error >> clCreateCommandQueue");

//===========================//  
//  Load OpenCL Kernel File  //
//===========================//

  char *cl_filename = "trace.cl";
  char *kername = "trace_kernel";
  char *source_str;
  size_t source_size;
    
  LoadKernel(cl_filename, &source_str, &source_size);    

// Create & Buid Program and Kernel //
  cl_program program;
  cl_kernel kernel;

  SetProgKernel(&program, &kernel, context, source_str, source_size, 
              Devices, dev_sel, kername);

//  Getptxcode(program);

// MakeParticle Setting
  ANGLE pAng, pAzimuth;
  POSITION sig0pos[3];
  DIRECTION sig0dir[3];
  particle *p;
  p = (particle *)malloc(sizeof(particle));
  
  pAng = XYANGLE * M_PI / AS_REAL(180.0);
  pAzimuth = AZIMUTH * M_PI / AS_REAL(180.0);
//  R[CX] = SIN(pAzimuth) * COS(pAng);
//  R[CY] = SIN(pAzimuth) * SIN(pAng);
//  R[CZ] = COS(pAzimuth);
//  for (i=0; i<3; i++) {
//      sig0pos[i] = ARENA->CEN[i] + R[i] * ARENA->RADIUS;
//      sig0dir[i] = -R[i];
//  }
    position_particle(pAzimuth, pAng, p);
    sig0pos[CX] = p->pos[CX];
    sig0pos[CY] = p->pos[CY];
    sig0pos[CZ] = p->pos[CZ];
    sig0dir[CX] = p->dir[CX];
    sig0dir[CY] = p->dir[CY];
    sig0dir[CZ] = p->dir[CZ];

//==============//  
// INPUT Buffer //
//==============//
//  unsigned int ran_seed = 90862346; 
#if L_TBL     
  cl_event IkernelEvent;
  size_t one = 1;
  size_t maxCU = 100;
#endif

  size_t gridsize = nx * ny * nz;
  int i;
  
  int ctrlflag = 0;  // 11 : initial move global table to local
                     // 66 & 67: run montecalo kernel
                     
  unsigned int *ran_seed;   
  
  ran_seed = (unsigned int *)malloc((int)ItermsPerBatch * sizeof(unsigned int));


  for (i = 0; i < ItermsPerBatch; i++)
  {
     ran_seed[i] = lrand48();
  }
  	  
  // control_flag   
  cl_mem control_flag = clCreateBuffer(context, 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
      sizeof(int), &ctrlflag, &ret);
  cl_errChk(ret,"Error in_seed");

      
  // rand seed 
  cl_mem in_seed = clCreateBuffer(context, 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
      ItermsPerBatch * sizeof(unsigned int), &ran_seed[0], &ret);
  cl_errChk(ret,"Error in_seed");
  
  int g_seed_size = (int)ItermsPerBatch * sizeof(unsigned int);
  
  // sigma0 pos 
  cl_mem cl_sig0pos = clCreateBuffer(context, 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
      sizeof(sig0pos), &sig0pos[0], &ret);
  cl_errChk(ret,"Error cl_sig0pos");
  
  int g_sig0pos_size = sizeof(sig0pos);

  ret = clEnqueueWriteBuffer(command_queue, cl_sig0pos, CL_TRUE, 0,
            sizeof(sig0pos), &sig0pos[0], 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_sig0pos");

  // sigma0 dir 
  cl_mem cl_sig0dir = clCreateBuffer(context, 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
      sizeof(sig0dir), &sig0dir[0], &ret);
  cl_errChk(ret,"Error cl_sig0dir");
  
  int g_sig0dir_size = sizeof(sig0dir);

  ret = clEnqueueWriteBuffer(command_queue, cl_sig0dir, CL_TRUE, 0,
            sizeof(sig0dir), &sig0dir[0], 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_sig0dir");
  
  // arena struct 
#if UNIFORM
  cl_mem cl_arena = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(arena), ARENA, &ret);
  cl_errChk(ret,"Error cl_arena");
  
  *g_arena_size = sizeof(arena);

  ret = clEnqueueWriteBuffer(command_queue, cl_arena, CL_TRUE, 0,
            sizeof(arena), 
            ARENA, 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_arena");
#else
  cl_mem cl_arena = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(arena) + (nx*ny*nz - 1)*sizeof(short), ARENA, &ret);
  cl_errChk(ret,"Error cl_arena");
  
  *g_arena_size = sizeof(arena);

  ret = clEnqueueWriteBuffer(command_queue, cl_arena, CL_TRUE, 0,
            sizeof(arena) + (nx*ny*nz - 1)*sizeof(short), 
            ARENA, 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_arena");
#endif
  
// ENRG  
  cl_mem cl_ENRG = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          NE * sizeof(TBL_REAL),
          ENRG, &ret);
  cl_errChk(ret,"Error cl_ENRG");
  
  int c_ENRG_size = NE * sizeof(TBL_REAL);
// write ENRG  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ENRG, CL_TRUE, 0,
          NE * sizeof(TBL_REAL),
          ENRG , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ENRG");

// PERC  
  cl_mem cl_PERC = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          NQ * sizeof(TBL_REAL),
          PERC, &ret);
  cl_errChk(ret,"Error cl_PERC");
  
  int c_PERC_size = NQ * sizeof(TBL_REAL);
// write PERC  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_PERC, CL_TRUE, 0,
          NQ * sizeof(TBL_REAL),
          PERC , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_PERC");

// FLUC  
  cl_mem cl_FLUC = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          NE*NQ*NM * sizeof(TBL_REAL),
          FLUC, &ret);
  cl_errChk(ret,"Error cl_FLUC");
  
  int c_FLUC_size = NE*NQ*NM * sizeof(TBL_REAL);
// write FLUC  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_FLUC, CL_TRUE, 0,
          NE*NQ*NM * sizeof(TBL_REAL),
          FLUC , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_FLUC");

// ProtonAdiposeTissue_Energy  
  cl_mem cl_ProtonAdiposeTissue_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonAdiposeTissue_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonAdiposeTissue_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonAdiposeTissue_Energy");
  
  int c_ProtonAdiposeTissue_Energy_size = sizeof(table) + (ProtonAdiposeTissue_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonAdiposeTissue_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonAdiposeTissue_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonAdiposeTissue_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonAdiposeTissue_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonAdiposeTissue_Energy");

// ProtonATissue_Energy  
  cl_mem cl_ProtonATissue_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonATissue_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonATissue_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonATissue_Energy");
  
  int c_ProtonATissue_Energy_size = sizeof(table) + (ProtonATissue_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonATissue_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonATissue_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonATissue_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonATissue_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonATissue_Energy");

// ProtonMuscleWithSucrose_Energy  
  cl_mem cl_ProtonMuscleWithSucrose_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonMuscleWithSucrose_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonMuscleWithSucrose_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonMuscleWithSucrose_Energy");
  
  int c_ProtonMuscleWithSucrose_Energy_size = sizeof(table) + (ProtonMuscleWithSucrose_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonMuscleWithSucrose_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonMuscleWithSucrose_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonMuscleWithSucrose_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonMuscleWithSucrose_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonMuscleWithSucrose_Energy");

// ProtonBBone_Energy  
  cl_mem cl_ProtonBBone_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonBBone_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonBBone_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonBBone_Energy");
  
  int c_ProtonBBone_Energy_size = sizeof(table) + (ProtonBBone_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonBBone_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonBBone_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonBBone_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonBBone_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonBBone_Energy");

// ProtonWater_Energy  
  cl_mem cl_ProtonWater_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonWater_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonWater_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonWater_Energy");
  
  *c_ProtonWater_Energy_size = sizeof(table) + (ProtonWater_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonWater_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonWater_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonWater_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonWater_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonWater_Energy");


// ProtonAir_Energy  
  cl_mem cl_ProtonAir_Energy = clCreateBuffer(context, 
          CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(table) + (ProtonAir_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonAir_Energy, &ret);
  cl_errChk(ret,"Error cl_ProtonAir_Energy");
  
  int c_ProtonAir_Energy_size = sizeof(table) + (ProtonAir_Energy->np - 1) * sizeof(TBL_REAL);
// write ProtonAir_Energy  <read only>
  ret = clEnqueueWriteBuffer(command_queue, cl_ProtonAir_Energy, CL_TRUE, 0,
          sizeof(table) + (ProtonAir_Energy->np - 1) * sizeof(TBL_REAL),
          ProtonAir_Energy , 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write cl_ProtonAir_Energy");
    
//===============//  
// OUTPUT Buffer //
//===============//
// Create Output memories //      
  // trace struct <dual buffer>
    cl_mem cl_traceA = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
            ItermsPerBatch * sizeof(trace),
            NULL, &ret);
    cl_errChk(ret,"Error >> cl_trace0");

    *g_traceA_size = ItermsPerBatch * sizeof(trace);
    
    cl_mem cl_traceB = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
            ItermsPerBatch * sizeof(trace),
            NULL, &ret);
    cl_errChk(ret,"Error >> cl_trace1");        
    
    int g_traceB_size = ItermsPerBatch * sizeof(trace);
    
    //=========//  
    // Mapping //
    //=========//
    
    // Mapping Kernel Parameters  //
      int a=0;
      ret = CL_SUCCESS;
      REAL CL_EMEAN, CL_ESTDV, CL_SIGMA, CL_XYANGLE, CL_AZIMUTH;
      CL_EMEAN   = (REAL) EMEAN;
      CL_ESTDV   = (REAL) ESTDV;
      CL_SIGMA   = (REAL) SIGMA;
      CL_XYANGLE = (REAL) XYANGLE;
      CL_AZIMUTH = (REAL) AZIMUTH;
        
//      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&control_flag);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&MAXSTEP);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&CL_EMEAN);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&CL_ESTDV);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&CL_SIGMA);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&CL_XYANGLE);
      ret |= clSetKernelArg(kernel, a++, sizeof(REAL), (void *)&CL_AZIMUTH);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&gridsize);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ItermsPerBatch);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_sig0pos);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_sig0dir);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&in_seed);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_traceA);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_traceB);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_arena);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&NE);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&NQ);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&NM);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ENRG);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_PERC);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_FLUC);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonWater_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonWater_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonAir_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonAir_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonAdiposeTissue_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonAdiposeTissue_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonATissue_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonATissue_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonMuscleWithSucrose_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonMuscleWithSucrose_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(int), (void *)&ProtonBBone_Energy->np);
      ret |= clSetKernelArg(kernel, a++, sizeof(table) + (ProtonBBone_Energy->np - 1) * sizeof(TBL_REAL), NULL);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonWater_Energy);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonAir_Energy);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonAdiposeTissue_Energy);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonATissue_Energy);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonMuscleWithSucrose_Energy);
      ret |= clSetKernelArg(kernel, a++, sizeof(cl_mem), (void *)&cl_ProtonBBone_Energy);
  
      cl_errChk(ret,"Error in Kernel Mapping");    

  int BatCnt;
  trace *cl_TraceOut;
  cl_TraceOut = (trace *)malloc(ItermsPerBatch * sizeof(trace));

  printf("Total Batches = %zd\n", Batches);

// loop zero => initializaion
  
#if L_TBL        
  ctrlflag = 11;

  ret = clEnqueueWriteBuffer(command_queue, control_flag, CL_TRUE, 0,
            sizeof(int), &ctrlflag, 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write control_flag");
  
      //===============================//  
      // Show Time for Initialization !!  //
      //===============================//
        
  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
      &maxCU, &one, 0, NULL, &IkernelEvent);
  cl_errChk(ret,"Error clEnqueueNDRangeKernel");
  clFinish(command_queue);
  ret = clGetEventProfilingInfo(IkernelEvent,CL_PROFILING_COMMAND_START,
          sizeof(cl_ulong),&eventStart,NULL);
  ret = clGetEventProfilingInfo(IkernelEvent,CL_PROFILING_COMMAND_END,
          sizeof(cl_ulong),&eventEnd,NULL);
  initime += (eventEnd-eventStart);
#endif          
//==============================================================//

  
// loop 2->N
  
  ret = clEnqueueWriteBuffer(command_queue, in_seed, CL_TRUE, 0,
            ItermsPerBatch * sizeof(unsigned int), &ran_seed[0], 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write in_seed");
  
  ctrlflag = 66;
  ret = clEnqueueWriteBuffer(command_queue, control_flag, CL_TRUE, 0,
            sizeof(int), &ctrlflag, 0, NULL, &wrRanEvent);
  cl_errChk(ret,"Error Write control_flag");   

  if (0 == ItermsPerWGroup)
  {
      ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
          &ItermsPerBatch, NULL, 0, NULL, &MkernelEvent);
  }
  else
  {
      ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
          &ItermsPerBatch, &ItermsPerWGroup, 0, NULL, &MkernelEvent);
  }

  cl_errChk(ret,"Error clEnqueueNDRangeKernel");
  clFinish(command_queue);
  ret = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_START,
          sizeof(cl_ulong),&eventStart,NULL);
  ret = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_END,
          sizeof(cl_ulong),&eventEnd,NULL);
  kerneltime += (eventEnd-eventStart);
  
#if EVENT_STAMP        
  printf("MkernelEvent eventStart = %lu\n", eventStart);
  printf("MkernelEvent eventEnd = %lu\n", eventEnd);
#endif          

  for (BatCnt = 1; BatCnt < Batches; BatCnt++)
  {
//      printf("BatCnt = %d\n", BatCnt);

      ctrlflag = ctrlflag ^ 1;   //66 <-> 67   
      
      // do collection
      if (ctrlflag == 67)
      {
          ret = clEnqueueReadBuffer(cmdQ_Read, cl_traceA, CL_TRUE , 0,
                  ItermsPerBatch * sizeof(trace),
                  cl_TraceOut, 0, NULL, &rdEvent_trace);
      }
      if (ctrlflag == 66)
      {
          ret = clEnqueueReadBuffer(cmdQ_Read, cl_traceB, CL_TRUE , 0,
                  ItermsPerBatch * sizeof(trace),
                  cl_TraceOut, 0, NULL, &rdEvent_trace);
      }
      cl_errChk(ret,"Error rdEvent_trace");

      ret = clEnqueueWriteBuffer(command_queue, control_flag, CL_TRUE, 0,
                sizeof(int), &ctrlflag, 0, NULL, &wrRanEvent);
      cl_errChk(ret,"Error Write control_flag");

      //===============================//  
      // Show Time for Monte Carlo !!  //
      //===============================//
      if (0 == ItermsPerWGroup)
      {
          ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
              &ItermsPerBatch, NULL, 0, NULL, &MkernelEvent);
      }
      else
      {
          ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
              &ItermsPerBatch, &ItermsPerWGroup, 0, NULL, &MkernelEvent);
      }
      cl_errChk(ret,"Error clEnqueueNDRangeKernel");

      
      clFinish(cmdQ_Read);
      for (i = 0; i < ItermsPerBatch; i++)
      {
	        collect(C, p, &cl_TraceOut[i]);
      }
      clFinish(command_queue);
      
      ret = clGetEventProfilingInfo(rdEvent_trace,CL_PROFILING_COMMAND_START,
              sizeof(cl_ulong),&eventStart,NULL);
      ret = clGetEventProfilingInfo(rdEvent_trace,CL_PROFILING_COMMAND_END,
              sizeof(cl_ulong),&eventEnd,NULL);
      rdbacktime += (eventEnd-eventStart);            
#if EVENT_STAMP 
  printf("rdEvent_trace eventStart = %lu\n", eventStart);
  printf("rdEvent_trace eventEnd = %lu\n", eventEnd);
#endif
      
      ret = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_START,
              sizeof(cl_ulong),&eventStart,NULL);
      ret = clGetEventProfilingInfo(MkernelEvent,CL_PROFILING_COMMAND_END,
              sizeof(cl_ulong),&eventEnd,NULL);
      kerneltime += (eventEnd-eventStart);
#if EVENT_STAMP 
  printf("MkernelEvent eventStart = %lu\n", eventStart);
  printf("MkernelEvent eventEnd = %lu\n", eventEnd);
#endif
                  
  }
 
  if (ctrlflag == 66)
  {
      ret = clEnqueueReadBuffer(cmdQ_Read, cl_traceA, CL_TRUE, 0,
              ItermsPerBatch * sizeof(trace),
              &cl_TraceOut[0], 0, NULL, &rdEvent_trace);
  }
  if (ctrlflag == 67)
  {
      ret = clEnqueueReadBuffer(cmdQ_Read, cl_traceB, CL_TRUE, 0,
              ItermsPerBatch * sizeof(trace),
              &cl_TraceOut[0], 0, NULL, &rdEvent_trace);
  }
  cl_errChk(ret,"Error rdEvent_trace");
      
  ret = clGetEventProfilingInfo(rdEvent_trace,CL_PROFILING_COMMAND_START,
          sizeof(cl_ulong),&eventStart,NULL);
  ret = clGetEventProfilingInfo(rdEvent_trace,CL_PROFILING_COMMAND_END,
          sizeof(cl_ulong),&eventEnd,NULL);
  rdbacktime += (eventEnd-eventStart);      
#if EVENT_STAMP 
  printf("rdEvent_trace eventStart = %lu\n", eventStart);
  printf("rdEvent_trace eventEnd = %lu\n", eventEnd);
#endif
  
  cl_errChk(ret,"Error rdEvent_trace");
  
  clFinish(cmdQ_Read);
  
  for (i = 0; i < ItermsPerBatch; i++)
  {
	    collect(C, p, &cl_TraceOut[i]);
  }

  printf("Initialization time = %f sec\n", (float)initime / 1.e9);
  printf("Monte Carlo time = %f sec\n", (float)kerneltime / 1.e9);
  printf("Read back time = %f sec\n", (float)rdbacktime / 1.e9);
  totalTime = initime+kerneltime+rdbacktime;
  
//  printf("Total kernel time = %f\n\n", (float)totalTime / 1.e9);

//==============//  
// Free buffer  //
//==============//

ret |= clFlush(command_queue);
ret |= clFinish(command_queue);
ret |= clReleaseKernel(kernel);
ret |= clReleaseProgram(program);
ret |= clReleaseMemObject(control_flag);
ret |= clReleaseMemObject(cl_sig0pos);
ret |= clReleaseMemObject(cl_sig0dir);
ret |= clReleaseMemObject(in_seed);
ret |= clReleaseMemObject(cl_arena);
ret |= clReleaseMemObject(cl_traceA);
ret |= clReleaseMemObject(cl_traceB);
ret |= clReleaseMemObject(cl_ENRG);
ret |= clReleaseMemObject(cl_PERC);
ret |= clReleaseMemObject(cl_FLUC);
ret |= clReleaseMemObject(cl_ProtonWater_Energy);
ret |= clReleaseMemObject(cl_ProtonAir_Energy);
ret |= clReleaseMemObject(cl_ProtonAdiposeTissue_Energy);
ret |= clReleaseMemObject(cl_ProtonATissue_Energy);
ret |= clReleaseMemObject(cl_ProtonMuscleWithSucrose_Energy);
ret |= clReleaseMemObject(cl_ProtonBBone_Energy);
ret |= clReleaseCommandQueue(command_queue);
ret |= clReleaseContext(context);
cl_errChk(ret,"Error in finish");    

*g_mem_size = g_seed_size + g_sig0pos_size + g_sig0dir_size + *g_arena_size
             + *g_traceA_size + g_traceB_size;

#if L_TBL        
*c_mem_size = c_ENRG_size + c_PERC_size + c_FLUC_size;
*l_mem_size = c_ProtonAdiposeTissue_Energy_size +            
                 c_ProtonATissue_Energy_size +
                 c_ProtonMuscleWithSucrose_Energy_size +
                 c_ProtonBBone_Energy_size +
                 *c_ProtonWater_Energy_size +
                 c_ProtonAir_Energy_size;
#else             
*c_mem_size = c_ENRG_size + c_PERC_size + c_FLUC_size + 
                 c_ProtonAdiposeTissue_Energy_size +
                 c_ProtonATissue_Energy_size +
                 c_ProtonMuscleWithSucrose_Energy_size +
                 c_ProtonBBone_Energy_size +
                 *c_ProtonWater_Energy_size +
                 c_ProtonAir_Energy_size;
*l_mem_size = 0;             
#endif          
}
