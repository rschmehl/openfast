// routines in KFAST_Lib.dll

#ifdef __cplusplus
#define EXTERNAL_ROUTINE extern "C"
#else
#define EXTERNAL_ROUTINE extern
#endif

EXTERNAL_ROUTINE void KFAST_Init(double *dt, int *numFlaps, int *numPylons, int *numComp, int numCompNds[], int modFlags[], const char KAD_FileName[], const char IfW_FileName[], const char MD_FileName[], const char KFC_FileName[],
    const char outFileRoot[], double *gravity, double windPt[], double FusODCM[], int *numRtrPtsElem, double rtrPts[], int *numRefPtElem, double refPts[], 
    int *numNodePtElem, double nodePts[], int *numDCMElem, double nodeDCMs[], 
    int *nFusOuts, int FusOutNd[], int *nSWnOuts, int SWnOutNd[], int *nPWnOuts, int PWnOutNd[], int *nVSOuts, int VSOutNd[], int *nSHSOuts, int SHSOutNd[], int *nPHSOuts, int PHSOutNd[], int *nPylOuts, int PylOutNd_c[], int *numOutChan, char* chanList[], int *errStat, char errMsg[]);
    // int *nFusOuts, int FusOutNd[], int *nSWnOuts, int SWnOutNd[], int *nPWnOuts, int PWnOutNd[], int *numOutChan, char chanList[][11], int *errStat, char errMsg[]);
EXTERNAL_ROUTINE void KFAST_AssRes(double *t, int* isInitialTime, int *numRtSpdRtrElem, double RtSpd_PyRtr[], double WindPt[], double FusO_prev[], double FusO[], double FusODCM_prev[], double FusODCM[], double FusOv_prev[],
                                    double FusOomegas_prev[], double FusOacc_prev[], int *numNodePtElem, double nodePts[],
                                    int *numNodeVelElem, double nodeVels[], int *numNodeOmegaElem, double nodeOmegas[], int *numNodeAccElem, double nodeAccs[],
                                    int *numDCMElem, double nodeDCMs[], int *numRtrPtsElem, double rtrPts[], double rtrVels[], double rtrDCMs[], 
                                    int *numNodeLoadsElem, double nodeLoads[], int *numRtrLoadsElem, double rtrLoads[], int *errStat, char errMsg[]);
EXTERNAL_ROUTINE void KFAST_AfterPredict(int *errStat, char errMsg[]);
EXTERNAL_ROUTINE void KFAST_Output(double *t, int *numGaussPtLoadsElem, double gaussPtLoads[], int *errStat, char errMsg[]);
EXTERNAL_ROUTINE void KFAST_End(int *errStat, char errMsg[]);

// some constants (keep these synced with values in FAST's fortran code)
#define INTERFACE_STRING_LENGTH 1025

#define ErrID_None 0 
#define ErrID_Info 1 
#define ErrID_Warn 2 
#define ErrID_Severe 3 
#define ErrID_Fatal 4 

static int AbortErrLev = ErrID_Fatal;      // abort error level; compare with NWTC Library

#define CHANNEL_LENGTH 10  
