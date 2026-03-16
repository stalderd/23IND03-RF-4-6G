# A2.3.2

METAS, CMI, NPL, PTB, TUBITAK and VSL will apply the method to extender-based VNA systems up to 220 GHz and implement the necessary software to automate the process. At least two out of the four waveguide bands 50 GHz   75 GHz, 75 GHz – 110 GHz, 110 GHz – 170 GHz and 170 GHz – 220 GHz will be evaluated. VDI will provide technical advice on the use of the extenders for this purpose.

## VNA Linearity Simulation Quadratic Model

The MATLAB script called [AnalyzeData_02.m](VNA_Linearity_Simulation_Quadratic_Model/AnalyzeData_02.m) analyze the simulated data sets `LinTestDataSet*.mat` and generates the following report, see [Simulation_opt.pdf](VNA_Linearity_Simulation_Quadratic_Model/Simulation_opt.pdf).

The purpose of the above script using simulated test data is to test if the quadratic model algorithm is working correctly.

## VNA Linearity Simulation Uncompressed Receiver Calibration Method

The MATLAB script called [AnalyzeData_UncompRecCalMet.m](VNA_Linearity_Simulation_Uncompressed_Receiver_Calibration_Method/AnalyzeData_UncompRecCalMet.m) analyzes the simulated data sets `LinTestDataSet*.mat` and generates the following report, see [Simulation_StepAtt_UncompresseReceiverCalMethod.pdf](VNA_Linearity_Simulation_Uncompressed_Receiver_Calibration_Method/Simulation_StepAtt_UncompresseReceiverCalMethod.pdf).

The purpose of the above script using simulated test data is to test if the uncompressed receiver calibration algorithm is working correctly.

## VNA Linearity Quadratic Model

- The MATLAB script called [AnalyzeData_opt_P1_01.m](VNA_Linearity_Quadratic_Model/AnalyzeData_opt_P1_01.m) analyzes the measured data sets for a1_p1 and b2_p1 and generates a report.
- The MATLAB script called [AnalyzeData_opt_P2_01.m](VNA_Linearity_Quadratic_Model/AnalyzeData_opt_P2_01.m) analyzes the measured data sets for a2_p2 and b1_p2 and generates a report.
- The MATLAB script called [E5061B_AnalyzeData_opt_P1_01.m](VNA_Linearity_Quadratic_Model/E5061B_AnalyzeData_opt_P1_01.m) analyzes the measured data sets for a1_p1 and b2_p1 and generates a report for a E5061B.
- The MATLAB script called [E5061B_AnalyzeData_opt_P2_01.m](VNA_Linearity_Quadratic_Model/E5061B_AnalyzeData_opt_P2_01.m) analyzes the measured data sets for a2_p2 and b1_p2 and generates a report for a E5061B.

## VNA Linearity Uncompressed Receiver Calibration Method

- The MATLAB script called [AnalyzeData_URCM_P1_01.m](VNA_Linearity_Uncompressed_Receiver_Calibration_Method/AnalyzeData_URCM_P1_01.m) analyzes the measured data sets and generates a report.
- The MATLAB script called [E5061B_AnalyzeData_URCM_P1_01.m](VNA_Linearity_Uncompressed_Receiver_Calibration_Method/E5061B_AnalyzeData_URCM_P1_01.m) analyzes the measured data sets and generates a report for a E5061B.
