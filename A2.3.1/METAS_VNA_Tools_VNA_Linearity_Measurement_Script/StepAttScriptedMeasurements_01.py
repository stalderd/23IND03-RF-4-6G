# Peter Huerlimann METAS - 22.03.2019
# Michael Wollensack METAS - 11.03.2022
# Daniel Stalder METAS - 2026

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
vnaDevice = 'Agilent_PNA_N5227A_#1' # Database file name of the VNA used in VNA Tools without the file extension (.vnadev)
sweepPoints = 100 # Number of sweep points per frequency point (CW sweep mode)
dwellTime = 0 # VNA dwell time (e.g. 0 s)
frequencies = [50e6, 1e9, 10e9, 18e9, 30e9, 40e9, 50e9] # VNA frequencies to be measured
sourcePowerLevels = range(-30, -4, 1) # VNA source power levels to be measured (e.g. -30 dBm to -5 dBm in 1 dB steps)

iSwitch = s.OpenSwitch(r'Agilent11713A', r'GPIB0::8::INSTR') # Step attenuator VNA Tools switch driver and visa resource name
stepAttAttenuation = range(0, 70, 10) # Step attenuator values (only used for naming, e.g. 0 dB to 60 dB it 10 dB steps)
stepAttStates = ('000','100','010','001','101','011','111') # Step attenuator states to be measured (e.g. 84905M)

vnaPowSetDelay = 100 # Delay after VNA source power level changes in ms (e.g. 100 ms)
vnaFreqSetDelay = 300000 # Delay after VNA frequency changes and after changing the source power from max to min in ms (e.g. 300'000 ms = 5 min) -> can be important for the receiver stability to wait long enough
stepAttDelay = 500 # Delay after step attenuator switching in ms (e.g. 500 ms)

directory = 'Measurements_01' # Measurement directory name
name = 'StepAtt60dB(f-f)_SN123456_01' # Measurement file base name (+ stepAttAttenuation[i].ToString() + 'dB.vdatb')
journalPath = 'Journal_01.vnalog' # Measurement journal file name (measurement journal must already exist)


journal = s.LoadJournal(journalPath)
journal.AddUserCommentJournalItem('Start scripted measurements')

iVna = s.OpenVna(vnaDevice)
mode = VnaParameter.SParameterMatrixReferenceReceiver(Array[int]([1, 2]))
iVna.ParameterMatrix = mode
iVna.SweepMode = VnaSweepMode.CWTime
iVna.SweepPoints = sweepPoints
iVna.DwellTime = dwellTime

for frequency in frequencies:
    print 'Frequency: ' + frequency.ToString() + ' Hz'
    
    fdir = directory + '\\' + name + '\\' + frequency.ToString() + 'Hz'
    Directory.CreateDirectory(RootPath + '\\' + fdir)

    iVna.FrequencyCW = frequency
    iVna.Source1Power = sourcePowerLevels[0]
    iVna.Source2Power = sourcePowerLevels[0]
    iVna.OutputState = True
    Thread.Sleep(vnaFreqSetDelay)
    
    for i in range(0, len(stepAttAttenuation)):
        print 'Attenuation: ' + stepAttAttenuation[i].ToString() + ' dB'
        
        iVna.Source1Power = sourcePowerLevels[0]
        iVna.Source2Power = sourcePowerLevels[0]

        journal.AddVnaSettingsJournalItem(vnaDevice, iVna)
        
        iSwitch.SetState(stepAttStates[i])
        Thread.Sleep(stepAttDelay)
        
        for sourcePower in sourcePowerLevels:
            print 'Source Power: ' + sourcePower.ToString() + ' dBm'
            
            pdir = fdir + '\\' + sourcePower.ToString() + 'dBm'
            Directory.CreateDirectory(RootPath + '\\' + pdir)
            
            iVna.Source1Power = sourcePower
            iVna.Source2Power = sourcePower
            journal.AddVnaSettingsJournalItem(vnaDevice, iVna)
            
            Thread.Sleep(vnaPowSetDelay)
            
            filename = pdir + '\\' + name + '_' + stepAttAttenuation[i].ToString() + 'dB.vdatb'
            iVna.TriggerSingle()
            data = iVna.GetData(VnaFormat.RawData)
            Binary.SaveVnaData[Number](s.RootPath + '\\' + filename, data)
            journal.AddMeasurementJournalItem(vnaDevice, filename)

iVna.Source1Power = sourcePowerLevels[0]
iVna.Source2Power = sourcePowerLevels[0]
iSwitch.SetState(stepAttStates[0])

iSwitch.Close()
iVna.Close()

journal.AddUserCommentJournalItem('End scripted measurements')

s.SaveJournal(journalPath, journal)
