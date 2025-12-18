# Peter Huerlimann METAS - 22.03.2019
# Michael Wollensack METAS - 11.03.2022
# Daniel Stalder METAS - 27.03.2025

# METAS VNA Tools project must already exist with Jounal before this script can be started.
# VNA must already be set correctly for the extender setup.
# IFBW must already be set correctly (e.g. 100 Hz).
# Point AVG to 1.

import clr

clr.AddReference('System.Windows.Forms')
clr.AddReference('Metas.Vna.Tools')
clr.AddReference('Metas.Vna.Journal')
clr.AddReference('Metas.Vna.Data')
clr.AddReference('Metas.UncLib.Core')
clr.AddReference('Metas.UncLib.LinProp')
clr.AddReference('Metas.Instr.Driver')

from System import Array
from System.IO import Directory
from System.Threading import Thread
from System.Windows.Forms import MessageBox, MessageBoxButtons
from Metas.Vna.Tools import Script
from Metas.Vna.Data import VnaParameter
from Metas.Vna.Data.FileIO import Binary
from Metas.UncLib.Core import Number
from Metas.UncLib.LinProp import UncNumber
from Metas.Instr.Driver.Vna import VnaSweepMode, VnaFormat

# Settings
vnaDevice = 'Keysight_PNA_N5227B_#2_VDI_WR10_2' # Database file name of the VNA used in VNA Tools without the file extension (.vnadev)
sweepPoints = 40 # Number of sweep points per frequency point (CW sweep mode)
vnaSetDelay_ms = 0 # Delay after a frequency change in ms

frequencies_Hz = range(int(68e9), int(118e9), int(2e9)) # Frequencies to be measured: start (inclusive), stop (exculsive), step

sourceAttenuator_dB = [0, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] # Extender source attenuator values to be measured
sourceAttenuator_position_port_1 = [0, 0.68, 0.80, 0.96, 1.19, 1.35, 1.50, 1.64, 1.750, 1.850, 1.96, 2.06, 2.145, 2.245, 2.330, 2.405, 2.490, 2.555] # Corresponding extender 1 source attenuator micrometer positions
sourceAttenuator_position_port_2 = [0, 0.71, 0.84, 1.00, 1.23, 1.40, 1.55, 1.69, 1.805, 1.915, 2.02, 2.12, 2.215, 2.310, 2.395, 2.480, 2.565, 2.645] # Corresponding extender 2 source attenuator micrometer positions

DUT = ('Thru_01','6dB_236275','10dB_0006','20dB_236276','20dB_236276_10dB_0006', 'Thru_02') # DUTs to be measured (e.g. Thru connection, 6 dB Attenuator, 10 dB Attenuator, 20 dB Attenuator, 20 dB Attenuator + 10 dB Attenuator  

s = Script(RootPath)
directory = 'Measurements_01' # Measurement directory name
Directory.CreateDirectory(RootPath + '\\' + directory)
journalPath = 'Journal_01.vnalog' # Measurement journal file name (measurement journal must already exist)

journal = s.LoadJournal(journalPath)
journal.AddUserCommentJournalItem('Start scripted measurements')

iVna = s.OpenVna(vnaDevice)
# mode = VnaParameter.SParameterMatrixAndReferenceReceiver(Array[int]([1, 2]))

mode_P1 = Array.CreateInstance(VnaParameter, 3, 1)
mode_P1[0, 0] = VnaParameter.FromString('S1,1')
mode_P1[1, 0] = VnaParameter.FromString('S2,1')
mode_P1[2, 0] = VnaParameter.FromString('a1_p1')

mode_P2 = Array.CreateInstance(VnaParameter, 3, 1)
mode_P2[0, 0] = VnaParameter.FromString('S2,2')
mode_P2[1, 0] = VnaParameter.FromString('S1,2')
mode_P2[2, 0] = VnaParameter.FromString('a2_p2')

iVna.ParameterMatrix = mode_P1
iVna.SweepMode = VnaSweepMode.CWTime
iVna.SweepPoints = sweepPoints

for i_DUT in range(0, len(DUT)):
	iVna.ParameterMatrix = mode_P1
    MessageBox.Show('Connect: ' + DUT[i_DUT], 'New Connection', MessageBoxButtons.OK)
    print('DUT: ' + DUT[i_DUT])
    MessageBox.Show(
        'Extender 2 Micrometer Position: 1.69' + '\n\n'
		'Wait 20 minutes to allow thermal equilibrium to be reached.',
        'Set Extender 2 Attenuator to 10 dB',
        MessageBoxButtons.OK
    )
    i_Att = 0
    for attenuator_dB in sourceAttenuator_dB:
        MessageBox.Show(
            'Extender 1 Micrometer Position: ' + sourceAttenuator_position_port_1[i_Att].ToString("F3"),
            'Set Extender 1 Attenuator to ' + attenuator_dB.ToString("F1") +' dB',
            MessageBoxButtons.OK
        )
        print('Source Attenuator: ' + attenuator_dB.ToString() + ' dB')
        i_Att += 1
        Thread.Sleep(3000)

        for frequency in frequencies_Hz:
            iVna.FrequencyCW = frequency
            journal.AddVnaSettingsJournalItem(vnaDevice, iVna)
            Thread.Sleep(vnaSetDelay_ms)

            print('Frequency: ' + frequency.ToString() + ' Hz')

            fdir = directory + '\\' + frequency.ToString() + 'Hz'
            pdir = fdir + '\\P1_' + attenuator_dB.ToString() + 'dB'
            filename = pdir + '\\' + DUT[i_DUT] + '.vdatb'
            
            Directory.CreateDirectory(RootPath + '\\' + fdir)
            Directory.CreateDirectory(RootPath + '\\' + pdir)

            iVna.TriggerSingle()
            data = iVna.GetData(VnaFormat.RawData)
            iVna.TriggerCont()
            Binary.SaveVnaData[Number](s.RootPath + '\\' + filename, data)
            journal.AddMeasurementJournalItem(vnaDevice, filename)

        s.SaveJournal(journalPath, journal)

    iVna.ParameterMatrix = mode_P2
    MessageBox.Show(
        'Extender 1 Micrometer Position: 1.64',
        'Set Extender 1 Attenuator to 10 dB',
        MessageBoxButtons.OK
    )
    i_Att = 0
    for attenuator_dB in sourceAttenuator_dB:
        MessageBox.Show(
            'Extender 2 Micrometer Position: ' + sourceAttenuator_position_port_2[i_Att].ToString("F3"),
            'Set Extender 2 Attenuator to ' + attenuator_dB.ToString("F1") +' dB',
            MessageBoxButtons.OK
        )
        print('Source Attenuator: ' + attenuator_dB.ToString() + ' dB')
        i_Att += 1
        Thread.Sleep(3000)

        for frequency in frequencies_Hz:
            iVna.FrequencyCW = frequency
            journal.AddVnaSettingsJournalItem(vnaDevice, iVna)
            Thread.Sleep(vnaSetDelay_ms)

            print('Frequency: ' + frequency.ToString() + ' Hz')

            fdir = directory + '\\' + frequency.ToString() + 'Hz'
            pdir = fdir + '\\P2_' + attenuator_dB.ToString() + 'dB'
            filename = pdir + '\\' + DUT[i_DUT] + '.vdatb'
            
            Directory.CreateDirectory(RootPath + '\\' + fdir)
            Directory.CreateDirectory(RootPath + '\\' + pdir)

            iVna.TriggerSingle()
            data = iVna.GetData(VnaFormat.RawData)
            iVna.TriggerCont()
            Binary.SaveVnaData[Number](s.RootPath + '\\' + filename, data)
            journal.AddMeasurementJournalItem(vnaDevice, filename)

        s.SaveJournal(journalPath, journal)

iVna.Close()
journal.AddUserCommentJournalItem('End scripted measurements')
s.SaveJournal(journalPath, journal)

MessageBox.Show('Reload the Journal in the Measurement Journal Tab!', 'Information', MessageBoxButtons.OK)
