# Peter Huerlimann METAS - 22.03.2019
# Michael Wollensack METAS - 11.03.2022

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
vnaDevice = 'Agilent_PNA_N5227A_#1'
sweepPoints = 100
dwellTime = 0
frequencies = [10e6, 20e6, 50e6, 100e6, 200e6, 500e6, 1e9, 2e9, 5e9, 10e9, 15e9, 20e9, 25e9, 30e9, 35e9, 40e9, 45e9, 50e9]
sourcePowerLevels = range(-30, 1, 1) # Source unleveled above +10 dBm
vnaSetDelay = 100

iSwitch = s.OpenSwitch(r'Agilent11713A', r'visa://l-217-01-11.ad.metas/GPIB0::8::INSTR')
stepAttAttenuation = range(0, 70, 5)
stepAttStates = ('0000','1000','0100','1100','0010','1010','0001','1001','0101','1101','0011','1011','0111','1111')
stepAttDelay = 600000

directory = 'Measurements_01'
name = 'StepAtt65dB(f-f)_TH61350416_01'
journalPath = 'Journal_01.vnalog'


journal = s.LoadJournal(journalPath)

journal.AddUserCommentJournalItem('Start scripted measurements')

iVna = s.OpenVna(vnaDevice)
mode = VnaParameter.SParameterMatrixReferenceReceiverAndSwitchTerms(Array[int]([1, 2]))
iVna.ParameterMatrix = mode
iVna.SweepMode = VnaSweepMode.CWTime
iVna.SweepPoints = sweepPoints
iVna.DwellTime = dwellTime

for frequency in frequencies:
    print 'Frequency: ' + frequency.ToString() + ' Hz'
    
    fdir = directory + '\\' + name + '\\' + frequency.ToString() + 'Hz'
    Directory.CreateDirectory(RootPath + '\\' + fdir)
    


    for i in range(0, len(stepAttAttenuation)):
        print 'Attenuation: ' + stepAttAttenuation[i].ToString() + ' dB'
        
        iVna.FrequencyCW = frequency
        iVna.Source1Power = sourcePowerLevels[0]
        iVna.Source2Power = sourcePowerLevels[0]
        iVna.OutputState = True
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
            
            Thread.Sleep(vnaSetDelay)
            
            filename = pdir + '\\' + name + '_' + stepAttAttenuation[i].ToString() + 'dB.vdatb'
            #journal.AddMeasurementJournalItem(vnaDevice, filename, iVna)
            iVna.TriggerSingle()
            data = iVna.GetData(VnaFormat.RawData)
            Binary.SaveVnaData[Number](s.RootPath + '\\' + filename, data)
            journal.AddMeasurementJournalItem(vnaDevice, filename)
            
            

iSwitch.SetState(stepAttStates[0])

iSwitch.Close()
iVna.Close()

journal.AddUserCommentJournalItem('End scripted measurements')

s.SaveJournal(journalPath, journal)
