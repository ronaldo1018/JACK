#include <stdio.h>
#include <stdlib.h>
#include "CL/cl.h"


#define MAX_SOURCE_SIZE (0x100000)
#define MAX_STEPS 5
#define MaxEnergy 201

#define CSV 1
#define PRT 1

extern cl_uint numPlatforms;
extern cl_uint numDevices;
extern cl_platform_id *platforms;
extern cl_device_id* Devices;

