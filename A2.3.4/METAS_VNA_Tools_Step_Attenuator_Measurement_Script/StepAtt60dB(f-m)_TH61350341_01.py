# Peter Morrissey METAS - 13.02.2026
# Daniel Stalder and Michael Wollensack METAS - 13.02.2026 16:36

# METAS VNA Tools project must already exist with Jounal before this script can be started.
# IFBW must already be set correctly (e.g. 10 Hz).
# Point AVG to 1.

import clr

clr.AddReference('System.Windows.Forms')
clr.AddReference('Metas.Vna.Tools')
clr.AddReference('Metas.Vna.Journal')
clr.AddReference('Metas.Vna.Data')
clr.AddReference('Metas.UncLib.Core')
clr.AddReference('Metas.Instr.Driver')

from System import Array
from System.IO import Directory
from System.Threading import Thread
from System.Windows.Forms import MessageBox
from Metas.Vna.Tools import Script
from Metas.Vna.Journal import Journal
from Metas.Vna.Data import VnaParameter
from Metas.Vna.Data.FileIO import Binary
from Metas.UncLib.Core import Number
from Metas.Instr.Driver.Vna import VnaSweepMode
from Metas.Instr.Driver.Vna import VnaFormat

s = Script(RootPath)

# Settings
vnaDevice = 'Keysight_PNA_N5225B_MY63058036' # Database file name of the VNA used in VNA Tools without the file extension (.vnadev)
dwellTime = 0 # VNA dwell time (e.g. 0 s)

iSwitch = s.OpenSwitch(r'Agilent11713A', r'GPIB0::8::INSTR') # Step attenuator VNA Tools switch driver and visa resource name

# measurementDescription, stepAttState, sourcePowerLevel, averageFactor
measurements = [ 
    ['01_0dB_-15dBm_AVG1', '000', -15, 1], 
    ['02_10dB_-15dBm_AVG1', '100', -15, 1], 
    ['03_20dB_-15dBm_AVG1', '010', -15, 1], 
    ['04_30dB_-15dBm_AVG1', '001', -15, 1],
    ['05_0dB_-15dBm_AVG1', '000', -15, 1], 
    ['06_30dB_-8dBm_AVG1', '001', -8, 1],
    ['07_40dB_-8dBm_AVG1', '101', -8, 1], 
    ['08_50dB_-8dBm_AVG10', '011', -8, 10], 
    ['09_60dB_-8dBm_AVG100','111', -8, 100],
    ['10_0dB_-15dBm_AVG1', '000', -15, 1]
]

stepAttDelay = 1000 # Delay after step attenuator switching (e.g. 1000 ms)

directory = 'Measurements_01\DUTs' # Measurement directory name
name = 'StepAtt60dB(f-m)_TH61350341_01' # Measurement file base name, change ending to _02, _03, _04 for new connection/measurements 
journalPath = 'Journal_01.vnalog' # Measurement journal file name (measurement journal must already exist)

Directory.CreateDirectory(RootPath + '\\' + directory + '\\' + name)

journal = s.LoadJournal(journalPath)

iVna = s.OpenVna(vnaDevice)

mode = VnaParameter.SParameterMatrixAndReferenceReceiver(Array[int]([1, 2]))
iVna.ParameterMatrix = mode
iVna.DwellTime = dwellTime
iVna.Source1Power = measurements[0][2]
iVna.Source2Power = measurements[0][2]
iVna.IFAverageFactor = measurements[0][3]
iVna.OutputState = True
journal.AddVnaSettingsJournalItem(vnaDevice, iVna)

for i in range(0, len(measurements)):
    
    iVna.Source1Power = measurements[i][2]
    iVna.Source2Power = measurements[i][2]
    iVna.IFAverageFactor = measurements[i][3]
    journal.AddVnaSettingsJournalItem(vnaDevice, iVna)
    
    iSwitch.SetState(measurements[i][1])
    Thread.Sleep(stepAttDelay)
    
    journal.AddNewSwitchStateJournalItem('Custom', int(1))
    journal.AddNewSwitchStateJournalItem('Custom', int(2))
    
    filename = directory + '\\' + name + '\\' + name + '_' + measurements[i][1].ToString() + '_' + measurements[i][0] + '.vdatb'
    iVna.TriggerSingle()
    data = iVna.GetData(VnaFormat.RawData)
    Binary.SaveVnaData[Number](s.RootPath + '\\' + filename, data)
    journal.AddMeasurementJournalItem(vnaDevice, filename)
    print(filename)

iSwitch.SetState(measurements[0][1])
    
journal.AddEndSwitchStateJournalItem()

iSwitch.Close()

iVna.Source1Power = measurements[0][2]
iVna.Source2Power = measurements[0][2]
iVna.IFAverageFactor = measurements[0][3]
journal.AddVnaSettingsJournalItem(vnaDevice, iVna)

iVna.Close()

s.SaveJournal(journalPath, journal)
