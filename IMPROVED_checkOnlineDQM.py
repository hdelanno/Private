from math import *
from dqmjson_online import *
from ROOT import TFile, gStyle, TCanvas, TH1, TH1F, TH2, TString, TObject
from optparse import OptionParser
import subprocess
import time
import datetime

import smtplib
import email
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.MIMEImage import MIMEImage
from email import Encoders


parser = OptionParser()
#parser.add_option("-p", "--process", dest="process", default="Run2011A-PromptReco-v4", help="Era and processing to consider")
#parser.add_option("-o", "--out", dest="out", default="CheckOnlineDQM.root", help="Output file") 
(options, args) = parser.parse_args()

ofile = open("ReportFromOnlineDQM.txt","w")
ofile_SMS = open("ReportFromOnlineDQM_SMS.txt","w")

#define lumi limits
highStatLumiCosmics = 30
highStatLumiCollisions = 15

highStatLumi = highStatLumiCollisions

# define problems
Detproblem = ''
StoNproblem = ''
DetproblemSMS = ''
StoNproblemSMS = ''
TkproblemSMS = ''
Tkproblem = ''
isProblem = 0
isProblemSMS = 0
RunNumber = 0
FillNumber = 0
BeamMode = ''
LumiSection = 0
RunKey = ''
RunType = ''
HltKey = ''
IsCollisionsRun = 0
BadChannels = 0
BadActiveChannels = 0
isFEDProblem = 0
isAllFEDsInError = 0
isChannelProblem = 0
isActiveChannelProblem = 0
isAPVeTimingProblem = 0
isDetFractionProblem = 0
isDetStoNProblem = 0
isTkProblem = 0
wasFEDProblem = 0
wasChannelProblem = 0
wasActiveChannelProblem = 0
wasAPVeTimingProblem = 0
wasDetFractionProblem = 0
wasDetStoNProblem = 0
wasTkProblem = 0
wasNoStripDataProblem = 0
FETimeDiffNames = ['FETimeDiffTECB' , 'FETimeDiffTECF' , 'FETimeDiffTIB' , 'FETimeDiffTOB' ]
FETimeDiffEntries = []
FETimeDiffMean = []
FETimeDiffRMS  = []
ProcessedEvents = 0
isActiveTOB      = False
isActiveTIDPlus  = False
isActiveTIDMinus = False
isActiveTIB      = False
isActiveTECPlus  = False
isActiveTECMinus = False
isNoStripDataProblem = 0
BadMajorityAddressesYMean = 0
BadMajorityAddressesFEDs = ''


def checksubdet(x, y, report, nlayer, value):
    global Detproblem
    global StoNproblem
    global DetproblemSMS
    global StoNproblemSMS
    if x == 1: det = 'TECMinus'
    if x == 2: det = 'TECPlus'
    if x == 3: det = 'TIB'
    if x == 4: det = 'TIDMinus'
    if x == 5: det = 'TIDPlus'
    if x == 6: det = 'TOB'    
    if report.GetName()=="detFractionReportMap":
        Detproblem += "Problem in Fraction of Good Modules in " + det + " Layer " + str(y) + " Value: %.4f \n" % float(value) 
        DetproblemSMS += "Problem in " + det + str(y) + ":%.3f\n" % float(value)
    else:
        if nlayer == 1:
            StoNproblem += "Problem in StoN in " + det + " Layer " + str(y) + " Value: %.4f \n" % float(value)
            StoNproblemSMS += "Problem S/N in " + det + str(y) + ":%.3f\n" % float(value)
        else:
            StoNproblem += "Problem in StoN in " + det + " Layer " + str(y) + " Value: %.4f \n" % float(value)
            StoNproblemSMS += "Problem S/N in " + det + str(y) + ":%.3f\n" % float(value)


def getmodulelist(moduledir):
    allmodules = []
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", moduledir, True)
    for something in data:
        allmodules.append(something)
        #        print something + " " + data[something]['value']
        allmodules.sort()
    return allmodules
            
def getstriphv(ls):
    if ( ls < highStatLumi ):
        return 0
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", "Info/EventInfo", True)
    summary = data['reportSummaryMap']['rootobj']
    hvon = 1
    for thisls in range(ls - highStatLumi , ls):
        hvon *= summary.GetBinContent(summary.FindBin( thisls , 19 )) * summary.GetBinContent(summary.FindBin( thisls , 20 )) * summary.GetBinContent(summary.FindBin( thisls , 21 )) * summary.GetBinContent(summary.FindBin( thisls , 22 ))
    return hvon

def getpixelhv(ls):
    if ( ls < highStatLumi ):
        return 0
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", "Info/EventInfo", True)
    summary = data['reportSummaryMap']['rootobj']
    hvon = 1
    for thisls in range(ls - highStatLumi , ls):
        hvon *= summary.GetBinContent(summary.FindBin( thisls , 16 )) * summary.GetBinContent(summary.FindBin( thisls , 17 ))
    return hvon

def getstripdata(ls):
    if (ls >= 2500):
        print "strip data: too high LS number"
        return
    #first find last filled bin
    lastfilledLS = 0
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", "SiStrip/MechanicalView", True)
    dataPlot = data['DataPresentInLS']['rootobj']
    for thisls in range(ls - 50 , ls):
        thiscontent = dataPlot.GetBinContent( thisls , 1 )
        if dataPlot.GetBinContent( thisls , 1 ) > 0:
            lastfilledLS = thisls

    #analyze activity
    global isActiveTOB     
    global isActiveTIDPlus 
    global isActiveTIDMinus
    global isActiveTIB     
    global isActiveTECPlus 
    global isActiveTECMinus
    #check if the detector is sending data (flag 1), or not (flag 0.01). To ignore missing LS in online DQM, also check empty bins (flag 0) - do not consider them as TRK not active!
    #TID+-
    for thisls in range(lastfilledLS - 1 , lastfilledLS):
        if dataPlot.GetBinContent( thisls , 4 ) == 1 or dataPlot.GetBinContent( thisls , 4 ) == 0:
            isActiveTIDMinus = True
        if dataPlot.GetBinContent( thisls , 5 ) == 1 or dataPlot.GetBinContent( thisls , 5 ) == 0:
            isActiveTIDPlus = True
    #other partitions            
    for thisls in range(lastfilledLS - 1 , lastfilledLS):
        if dataPlot.GetBinContent( thisls , 1 ) == 1 or dataPlot.GetBinContent( thisls , 1 ) == 0:
            isActiveTECMinus = True
        if dataPlot.GetBinContent( thisls , 2 ) == 1 or dataPlot.GetBinContent( thisls , 2 ) == 0:
            isActiveTECPlus = True
        if dataPlot.GetBinContent( thisls , 3 ) == 1 or dataPlot.GetBinContent( thisls , 3 ) == 0:
            isActiveTIB = True
        if dataPlot.GetBinContent( thisls , 6 ) == 1 or dataPlot.GetBinContent( thisls , 6 ) == 0:
            isActiveTOB = True
        
def getmodulelistdict(moduledir):
    allmodules = []
    allmodulevalues = []
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", moduledir, True)
    for something in data:
        allmodules.append(something)
        reason = int(data[something]['value'])
        saved_reason = ''
        if ( (reason >> 0) & 0x1 ): saved_reason += "| FED bad channel "
        if ( (reason >> 1) & 0x1 ): saved_reason += "| # of Digi "
        if ( (reason >> 2) & 0x1 ): saved_reason += "| # of Clusters "
        if ( (reason >> 3) & 0x1 ): saved_reason += "| Excluded FED channel "
        if ( (reason >> 4) & 0x1 ): saved_reason += "| DCS error "
        allmodulevalues.append(saved_reason)
#        print something + " " + saved_reason + " " + str(reason)
    modules = dict( zip ( allmodules , allmodulevalues ) )
    return modules
            
def getmodulereasonsdict(moduledir):
    allmodules = []
    allmodulereasons = []
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", moduledir, True)
    for something in data:
        allmodules.append(something)
        reason = int(data[something]['value'])
        allmodulereasons.append(reason)
    modules = dict( zip ( allmodules , allmodulereasons ) )
    return modules
            

def isOnlinePublishing():
    status = True
    data = dqm_get_json(serverurl, 0 , "/Online/ALL", "/SiStrip/EventInfo", True)
    if str(data) == "{}":
        status = False
    return status

def getTrivials():
    global RunNumber
    global LumiSection
    global RunKey
    global RunType
    global HltKey
    global IsCollisionsRun
    global BeamMode
    global FillNumber
    global ProcessedEvents
    global BadMajorityAddressesYMean
    global BadMajorityAddressesFEDs

    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", '/SiStrip/ReadoutView/FE/VsId', True)

    bmAddrs = data['BadMajorityAddresses']['rootobj']
    #    BadMajorityAddressesYMean = bmAddrs.GetMean(2)
    sum = 0
    nactivebins = 0
    for xbin in range(1, bmAddrs.GetNbinsX()):
        if bmAddrs.GetBinContent(xbin) > 0:
            sum += bmAddrs.GetBinContent(xbin)
            nactivebins += 1
            BadMajorityAddressesFEDs += 'FED' + str(int(bmAddrs.GetBinCenter(xbin))) + ' '
    #    BadMajorityAddressesYMean = sum / bmAddrs.GetNbinsX()
    if nactivebins > 0:
        BadMajorityAddressesYMean = sum / float(nactivebins)
    print "BadMajorityAddressesYMean = ", BadMajorityAddressesYMean
    if nactivebins > 5:
        BadMajorityAddressesFEDs = '0'
                            
    data = dqm_get_json(serverurl, 0 , "/Online/ALL", '/SiStrip/EventInfo', True)
    RunNumber = int(data['iRun']['value'])
    LumiSection = int(data['iLumiSection']['value'])
    ProcessedEvents = int(data['processedEvents']['value'])

    data = dqm_get_json(serverurl, 0 , "/Online/ALL", '/Info/ProvInfo', True)
    RunKey = str(data['Run Type']['value'])
    IsCollisionsRun = int(data['isCollisionsRun']['value'])
    HltKey = str(data['hltKey']['value'])
    if ( IsCollisionsRun == 1 ):
        RunType = 'Collisions'
    else:
        if ( RunKey.lower() == 'cosmic_run' ):
            RunType = 'Cosmics'
        else:
            RunType = 'Commissioning'
    
    #get beam mode
    data = dqm_get_json(serverurl, 0 , "/Online/ALL", '/Info/LhcInfo', True)
    beamsetupHisto = data['beamMode']['rootobj']
    BeamMode = beamsetupHisto.GetYaxis().GetBinLabel( int(beamsetupHisto.GetBinContent(LumiSection-1)) )
    if ( BeamMode.lower() == 'stable' ):
        BeamMode = 'stable beams'
    BeamMode = BeamMode.upper()
    lhcfillHisto = data['lhcFill']['rootobj']
    FillNumber = int(lhcfillHisto.GetBinContent ( LumiSection-1 ))

def getnFEDErrors():    
    fld = '/SiStrip/ReadoutView/FED'
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", fld, True)
    return float(data['nFEDErrors']['rootobj'].GetMean())

def getnBadChannels():
    global BadChannels
    global BadActiveChannels
    fld = '/SiStrip/ReadoutView/Fiber'
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", fld, True)
    BadChannels = float(data['nBadChannelStatusBits']['rootobj'].GetMean())
    BadActiveChannels = float(data['nBadActiveChannelStatusBits']['rootobj'].GetMean())
    

def getFETimeDiff(plot):
    fld = '/SiStrip/ReadoutView/FE/APVe'
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", fld, True)
    info = []
    info.append ( data[ str(plot) ]['rootobj'].GetMean() )
    info.append ( data[ str(plot) ]['rootobj'].GetRMS() )
    info.append ( data[ str(plot) ]['rootobj'].GetEntries() )
    return info

thistime = time.time()
localtime = time.localtime(thistime)
whattime = str(localtime[3])+":"+str(localtime[4])+":"+str(localtime[5])
#print whattime

#send heartbeat sms if in the time slot
sendHB_SMS = 0
if ( localtime[3] == 12 and localtime[4] > 29 and localtime[4] < 34 ):
    print "Time to send heartbeat SMS"
    #sendHB_SMS = 1

ofilehbsms = open("sendheartbeatsms.txt","w")
ofilehbsms.write(str(sendHB_SMS)+"\n")
ofilehbsms.close()
    
if isOnlinePublishing() == False:
    print "No online publishing."
    ofilee = open("sendemail.txt","w")
    ofilee.write(str(0)+"\n")
    ofilee.close()
    ofilesms = open("sendsms.txt","w")
    ofilesms.write(str(0)+"\n")
    ofilesms.close()
    sys.exit(0)

#first get runnumber, lumisection, etc
getTrivials()

if ( RunType.lower() != 'collisions' ):
    highStatLumi = 30

stripHv = getstriphv(LumiSection - 1)
pixelHv = getpixelhv(LumiSection - 1)

print "Run ", RunNumber, " LS ", LumiSection
print "Fill Number ", FillNumber
print "Strip HV ON for the past ", highStatLumi, "LS : ", stripHv
print "Pixel HV ON for the past ", highStatLumi, "LS : ", pixelHv
print "RunKey ", RunKey
print "IsCollisionsRun ", IsCollisionsRun
print "RunType ", RunType
print "HLTKey ", HltKey
print "Beam Mode ", BeamMode

#do not reprocess the same run again (for now)
#if ( RunNumber == old_run ):
#    sys.exit(0)

ofile_reference_previous = open("reference.txt","r")
old_run  = ofile_reference_previous.readline()
old_fill  = ofile_reference_previous.readline()
old_time = ofile_reference_previous.readline()
if ( int(RunNumber) == int(old_run) ):
    wasFEDProblem           = int(ofile_reference_previous.readline())
    wasChannelProblem       = int(ofile_reference_previous.readline())
    wasActiveChannelProblem = int(ofile_reference_previous.readline())
    wasAPVeTimingProblem    = int(ofile_reference_previous.readline())
    wasDetFractionProblem   = int(ofile_reference_previous.readline())
    wasDetStoNProblem       = int(ofile_reference_previous.readline())
    wasTkProblem            = int(ofile_reference_previous.readline())
    wasNoStripDataProblem   = int(ofile_reference_previous.readline())

#DEBUG
print "read wasNoStripDataProblem: ",wasNoStripDataProblem

#old_localtime = time.localtime(float(str(old_time)))

#read actual test limits
limits_filename = "limits_pprun.txt"
#RunKey = "cosmic_run"
if ( RunKey.lower() == 'cosmic_run' ):
    limits_filename = "limits_cosmicrun.txt"
limitfile = open( str(limits_filename) , "r" )

fed_limits = json.loads( limitfile.readline() )
apve_limits = json.loads( limitfile.readline() )
gooddet_limits = json.loads( limitfile.readline() )
sn_limits = json.loads( limitfile.readline() )
track_limits = json.loads( limitfile.readline() )
newbadmod_limits = json.loads( limitfile.readline() )
limitfile.close()

# Checking LS
ofile.write( "#####   Report from Online DQM for run "+ str(RunNumber) +"   ##### \n")
ofile.write("\n")
ofile.write( "BEAM MODE:   "+ str(BeamMode) + "\n")
ofile.write( "HLT KEY:     "+ str(HltKey) + "\n")
ofile.write( "RUN TYPE:    "+ str(RunType) + "\n")
ofile.write( "FILL NUMBER: "+ str(FillNumber) + "\n")
#ofile.write( "\nAnalysis done at LS "+ str(LumiSection) + " at " + whattime + "\n")
ofile.write( "\nAnalysis done at LS "+ str(LumiSection) + "\n")
ofile.write("\n")
ofile.write("-------------------------------------------- \n")

ofile_SMS.write( "R" + str(RunNumber) + ",LS" + str(LumiSection) + ":\n" )

getstripdata(LumiSection)

#DEBUG
print "conditions: ", ( isActiveTIDPlus == False or isActiveTIDMinus == False )
#rely only on TID now (low event rate in other partitions)
#if ( isActiveTOB == False or isActiveTIDPlus == False or isActiveTIDMinus == False or isActiveTIB == False or isActiveTECPlus == False or isActiveTECMinus == False ):
if ( wasNoStripDataProblem < 1 and LumiSection < 2500 and ( isActiveTIDPlus == False or isActiveTIDMinus == False or isActiveTOB == False or isActiveTIB == False or isActiveTECPlus == False or isActiveTECMinus == False )):
    #    detAffected = 'NO DIGIS for 10 LS in'
    detAffected = 'NO DIGIS in'
    if isActiveTOB == False:
        detAffected += " TOB"
    if isActiveTIDPlus == False:
        detAffected += " TID+"
    if isActiveTIDMinus == False:
        detAffected += " TID-"
    if isActiveTIB == False:
        detAffected += " TIB"
    if isActiveTECPlus == False:
        detAffected += " TEC+"
    if isActiveTECMinus == False:
        detAffected += " TEC-"
    if ( isActiveTOB == False and isActiveTIDPlus == False and isActiveTIDMinus == False and isActiveTIB == False and isActiveTECPlus == False and isActiveTECMinus == False ):
        detAffected += ": is HV on?"
    ofile.write(str(detAffected)+"\n")                
    ofile_SMS.write(str(detAffected)+"\n")
    isProblem = 1
    isProblemSMS = 1
    isNoStripDataProblem = 1

#DEBUG
print "isNoStripDataProblem: ", isNoStripDataProblem

nFEDErr = getnFEDErrors()
getnBadChannels()

if ( nFEDErr > 350 ):
    isAllFEDsInError = 1

if ( wasFEDProblem == 2 or ( wasFEDProblem == 1 and isAllFEDsInError == 0 ) ):
    nFEDErr = 0
    isAllFEDsInError = 0
if ( wasChannelProblem == 1 ):
    BadChannels = 0
if ( wasActiveChannelProblem == 1 ):
    BadActiveChannels = 0

if ((nFEDErr > int(fed_limits['federrors']) ) or (BadChannels > int(fed_limits['badchannels']) ) or (BadActiveChannels > int(fed_limits['badactivechannels'])) ):
    ofile.write( "Check: FED errors \n")
    ofile.write("\n")
    if ( nFEDErr > int(fed_limits['federrors']) ):
        ofile.write("Large number of FED with errors: %.2f\n" % float(nFEDErr))                
        ofile_SMS.write( "FEDErrors:%.2f\n" % float(nFEDErr))
        isFEDProblem = 1
        if ( isAllFEDsInError == 1 ):
            isFEDProblem = 2            
    if ( BadChannels > int(fed_limits['badchannels']) ):
        ofile.write("Large number of Bad Channels: " + str(BadChannels)+"\n")                
        ofile_SMS.write( "BadChannels:" + str(BadChannels) + "\n" )
        isChannelProblem = 1
    if ( BadActiveChannels > int(fed_limits['badactivechannels']) ):
        ofile.write("Large number of Bad Active Channels: " + str(BadActiveChannels)+"\n")                
        ofile_SMS.write( "BadActiveChannels:" + str(BadActiveChannels) + "\n" )
        isActiveChannelProblem = 1
    ofile.write("\n")
    ofile.write("-------------------------------------------- \n")
    isProblem = 1
    isProblemSMS = 1

if ( wasAPVeTimingProblem == 0 ):
    for apvetiming in FETimeDiffNames:
        thisTiming = getFETimeDiff( apvetiming )
        if ( ProcessedEvents > 100 and BadMajorityAddressesYMean/ProcessedEvents > float(apve_limits['apvetiming']) and int(thisTiming[2] ) > 0 ):
            #            if ( BadMajorityAddressesYMean/ProcessedEvents > float(apve_limits['apvetiming']) ):
            #if ( BadMajorityAddressesYMean > 1 ):

            ofile.write( "APVe mismatch for " + str(apvetiming) + ": mean = %.2f RMS = %.2f\n" % (float(thisTiming[0]) , float(thisTiming[1]) ))
            #                ofile_SMS.write( "APVeTiming:" + str(apvetiming) + "\n" )
            isAPVeTimingProblem = 1
            isProblem = 1
            #                isProblemSMS = 1
    if isAPVeTimingProblem == 1:
        ofile.write( " - observed in %.4f%% of processed events (~%i events)\n" % ( float(BadMajorityAddressesYMean/ProcessedEvents*100.), int(round(BadMajorityAddressesYMean)) ))
        if BadMajorityAddressesFEDs != '0':
            ofile.write( " - FEDs with bad majority APV addresses: " + BadMajorityAddressesFEDs + "\n" ) 
        ofile.write( "-------------------------------------------- \n" )
    
if ( stripHv == 1 and wasDetFractionProblem == 0 and RunType.lower() != 'commissioning' ):
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", '/SiStrip/MechanicalView', True)
    report1 = data['detFractionReportMap']['rootobj']
    nbadlayer1 = 0
    nbinX = report1.GetNbinsX()+1
    nbinY = report1.GetNbinsY()+1
    for xbin in range(1, nbinX):
        for ybin in range(1,nbinY):
            bin = report1.GetBin(xbin,ybin)
            dqm = report1.GetBinContent(bin)
            if dqm < -0.5: continue
            isbad = 0

            if ( xbin == 1 and ybin == 1 and dqm < float(gooddet_limits['tecb1']) ): isbad = 1
            if ( xbin == 1 and ybin == 2 and dqm < float(gooddet_limits['tecb2']) ): isbad = 1
            if ( xbin == 1 and ybin == 3 and dqm < float(gooddet_limits['tecb3']) ): isbad = 1
            if ( xbin == 1 and ybin == 4 and dqm < float(gooddet_limits['tecb4']) ): isbad = 1
            if ( xbin == 1 and ybin == 5 and dqm < float(gooddet_limits['tecb5']) ): isbad = 1
            if ( xbin == 1 and ybin == 6 and dqm < float(gooddet_limits['tecb6']) ): isbad = 1
            if ( xbin == 1 and ybin == 7 and dqm < float(gooddet_limits['tecb7']) ): isbad = 1
            if ( xbin == 1 and ybin == 8 and dqm < float(gooddet_limits['tecb8']) ): isbad = 1
            if ( xbin == 1 and ybin == 9 and dqm < float(gooddet_limits['tecb9']) ): isbad = 1
            if ( xbin == 2 and ybin == 1 and dqm < float(gooddet_limits['tecf1']) ): isbad = 1
            if ( xbin == 2 and ybin == 2 and dqm < float(gooddet_limits['tecf2']) ): isbad = 1
            if ( xbin == 2 and ybin == 3 and dqm < float(gooddet_limits['tecf3']) ): isbad = 1
            if ( xbin == 2 and ybin == 4 and dqm < float(gooddet_limits['tecf4']) ): isbad = 1
            if ( xbin == 2 and ybin == 5 and dqm < float(gooddet_limits['tecf5']) ): isbad = 1
            if ( xbin == 2 and ybin == 6 and dqm < float(gooddet_limits['tecf6']) ): isbad = 1
            if ( xbin == 2 and ybin == 7 and dqm < float(gooddet_limits['tecf7']) ): isbad = 1
            if ( xbin == 2 and ybin == 8 and dqm < float(gooddet_limits['tecf8']) ): isbad = 1
            if ( xbin == 2 and ybin == 9 and dqm < float(gooddet_limits['tecf9']) ): isbad = 1
            if ( xbin == 3 and ybin == 1 and dqm < float(gooddet_limits['tib1']) ): isbad = 1
            if ( xbin == 3 and ybin == 2 and dqm < float(gooddet_limits['tib2']) ): isbad = 1
            if ( xbin == 3 and ybin == 3 and dqm < float(gooddet_limits['tib3']) ): isbad = 1
            if ( xbin == 3 and ybin == 4 and dqm < float(gooddet_limits['tib4']) ): isbad = 1
            if ( xbin == 4 and ybin == 1 and dqm < float(gooddet_limits['tidb1']) ): isbad = 1
            if ( xbin == 4 and ybin == 2 and dqm < float(gooddet_limits['tidb2']) ): isbad = 1
            if ( xbin == 4 and ybin == 3 and dqm < float(gooddet_limits['tidb3']) ): isbad = 1
            if ( xbin == 5 and ybin == 1 and dqm < float(gooddet_limits['tidf1']) ): isbad = 1
            if ( xbin == 5 and ybin == 2 and dqm < float(gooddet_limits['tidf2']) ): isbad = 1
            if ( xbin == 5 and ybin == 3 and dqm < float(gooddet_limits['tidf3']) ): isbad = 1
            if ( xbin == 6 and ybin == 1 and dqm < float(gooddet_limits['tob1']) ): isbad = 1
            if ( xbin == 6 and ybin == 2 and dqm < float(gooddet_limits['tob2']) ): isbad = 1
            if ( xbin == 6 and ybin == 3 and dqm < float(gooddet_limits['tob3']) ): isbad = 1
            if ( xbin == 6 and ybin == 4 and dqm < float(gooddet_limits['tob4']) ): isbad = 1
            if ( xbin == 6 and ybin == 5 and dqm < float(gooddet_limits['tob5']) ): isbad = 1
            if ( xbin == 6 and ybin == 6 and dqm < float(gooddet_limits['tob6']) ): isbad = 1

            if ( isbad == 1 ):
                nbadlayer1 = nbadlayer1 + 1            
                checksubdet(xbin, ybin, report1, nbadlayer1, dqm)
                isDetFractionProblem = 1

if ( stripHv == 1 and IsCollisionsRun == 1 and wasDetStoNProblem == 0 ):
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", '/SiStrip/MechanicalView', True)
    report2 = data['sToNReportMap']['rootobj']
    nbadlayer2 = 0
    nbinX = report2.GetNbinsX()+1
    nbinY = report2.GetNbinsY()+1
    for xbin in range(1, nbinX):
        for ybin in range(1,nbinY):
            bin = report2.GetBin(xbin,ybin)
            dqm = report2.GetBinContent(bin)
            if dqm < -0.5: continue
            isbad = 0
            if ( xbin == 1 and ybin == 1 and dqm < float(sn_limits['sntecb1']) ): isbad = 1
            if ( xbin == 1 and ybin == 2 and dqm < float(sn_limits['sntecb2']) ): isbad = 1
            if ( xbin == 1 and ybin == 3 and dqm < float(sn_limits['sntecb3']) ): isbad = 1
            if ( xbin == 1 and ybin == 4 and dqm < float(sn_limits['sntecb4']) ): isbad = 1
            if ( xbin == 1 and ybin == 5 and dqm < float(sn_limits['sntecb5']) ): isbad = 1
            if ( xbin == 1 and ybin == 6 and dqm < float(sn_limits['sntecb6']) ): isbad = 1
            if ( xbin == 1 and ybin == 7 and dqm < float(sn_limits['sntecb7']) ): isbad = 1
            if ( xbin == 1 and ybin == 8 and dqm < float(sn_limits['sntecb8']) ): isbad = 1
            if ( xbin == 1 and ybin == 9 and dqm < float(sn_limits['sntecb9']) ): isbad = 1
            if ( xbin == 2 and ybin == 1 and dqm < float(sn_limits['sntecf1']) ): isbad = 1
            if ( xbin == 2 and ybin == 2 and dqm < float(sn_limits['sntecf2']) ): isbad = 1
            if ( xbin == 2 and ybin == 3 and dqm < float(sn_limits['sntecf3']) ): isbad = 1
            if ( xbin == 2 and ybin == 4 and dqm < float(sn_limits['sntecf4']) ): isbad = 1
            if ( xbin == 2 and ybin == 5 and dqm < float(sn_limits['sntecf5']) ): isbad = 1
            if ( xbin == 2 and ybin == 6 and dqm < float(sn_limits['sntecf6']) ): isbad = 1
            if ( xbin == 2 and ybin == 7 and dqm < float(sn_limits['sntecf7']) ): isbad = 1
            if ( xbin == 2 and ybin == 8 and dqm < float(sn_limits['sntecf8']) ): isbad = 1
            if ( xbin == 2 and ybin == 9 and dqm < float(sn_limits['sntecf9']) ): isbad = 1
            if ( xbin == 3 and ybin == 1 and dqm < float(sn_limits['sntib1']) ): isbad = 1
            if ( xbin == 3 and ybin == 2 and dqm < float(sn_limits['sntib2']) ): isbad = 1
            if ( xbin == 3 and ybin == 3 and dqm < float(sn_limits['sntib3']) ): isbad = 1
            if ( xbin == 3 and ybin == 4 and dqm < float(sn_limits['sntib4']) ): isbad = 1
            if ( xbin == 4 and ybin == 1 and dqm < float(sn_limits['sntidb1']) ): isbad = 1
            if ( xbin == 4 and ybin == 2 and dqm < float(sn_limits['sntidb2']) ): isbad = 1
            if ( xbin == 4 and ybin == 3 and dqm < float(sn_limits['sntidb3']) ): isbad = 1
            if ( xbin == 5 and ybin == 1 and dqm < float(sn_limits['sntidf1']) ): isbad = 1
            if ( xbin == 5 and ybin == 2 and dqm < float(sn_limits['sntidf2']) ): isbad = 1
            if ( xbin == 5 and ybin == 3 and dqm < float(sn_limits['sntidf3']) ): isbad = 1
            if ( xbin == 6 and ybin == 1 and dqm < float(sn_limits['sntob1']) ): isbad = 1
            if ( xbin == 6 and ybin == 2 and dqm < float(sn_limits['sntob2']) ): isbad = 1
            if ( xbin == 6 and ybin == 3 and dqm < float(sn_limits['sntob3']) ): isbad = 1
            if ( xbin == 6 and ybin == 4 and dqm < float(sn_limits['sntob4']) ): isbad = 1
            if ( xbin == 6 and ybin == 5 and dqm < float(sn_limits['sntob5']) ): isbad = 1
            if ( xbin == 6 and ybin == 6 and dqm < float(sn_limits['sntob6']) ): isbad = 1
            if ( isbad == 1 ):
                isDetStoNProblem = 1
                checksubdet(xbin,ybin,report2,nbadlayer2, dqm)
                nbadlayer2 = nbadlayer2 + 1

if ( ( isDetFractionProblem == 1 ) or ( isDetStoNProblem == 1 ) ):
    isProblem = 1
    isProblemSMS = 1
    ofile.write( "Check: Report Summary map \n")
    ofile.write("\n")
    ofile.write(Detproblem)
    ofile.write(StoNproblem)
    ofile.write("-------------------------------------------- \n")
    ofile_SMS.write(DetproblemSMS)
    ofile_SMS.write(StoNproblemSMS)
    

if ( stripHv == 1 and pixelHv == 1 and wasTkProblem == 0 and RunType.lower() == 'collisions' ):
    fld = '/Tracking/EventInfo'
    data = dqm_get_json(serverurl, str(RunNumber), "/Online/ALL", fld, True)
    tracking = data['reportSummaryMap']['rootobj']
    nbinX = tracking.GetNbinsX()+1
    nbinY = tracking.GetNbinsY()+1
    for xbin in range(1, nbinX):
        for ybin in range(1, nbinY):
            tkvalue = tracking.GetBinContent(xbin,ybin)
            if ( tkvalue < -0.5 and ( xbin == 1 or xbin == 2 or xbin == 3 ) ): continue        

            if ( xbin == 1 and tkvalue < float(track_limits['trackchi2']) ):
                Tkproblem += "Problem in Track reconstruction (chi2): %.4f\n" % float(tkvalue)
                TkproblemSMS += "TrackReco (chi2):%.3f\n" % float(tkvalue)
                isTkProblem = 1
            if ( xbin == 2 and tkvalue < float(track_limits['trackrate']) ):
                Tkproblem += "Problem in Track reconstruction (rate): %.4f\n" % float(tkvalue)
                TkproblemSMS += "TrackReco (rate):%.3f\n" % float(tkvalue)
                isTkProblem = 1
            if ( xbin == 3 and tkvalue < float(track_limits['trackrechits']) ):
                Tkproblem += "Problem in Track reconstruction (rechits): %.4f\n" % float(tkvalue)
                TkproblemSMS += "TrackReco (rechits):%.3f\n" % float(tkvalue)
                isTkProblem = 1

if ( isTkProblem == 1 ):
    isProblem = 1
    isProblemSMS = 1
    ofile.write("Check: Tracking report summary map \n")
    ofile.write("\n")
    ofile.write(Tkproblem)
    ofile.write("-------------------------------------------- \n")
    ofile_SMS.write(TkproblemSMS)


modulesTECMinus = getmodulelist('/SiStrip/MechanicalView/TEC/side_1/BadModuleList')
modulesTECPlus = getmodulelist('/SiStrip/MechanicalView/TEC/side_2/BadModuleList')
modulesTIB = getmodulelist('/SiStrip/MechanicalView/TIB/BadModuleList')
modulesTIDMinus = getmodulelist('/SiStrip/MechanicalView/TID/side_1/BadModuleList')
modulesTIDPlus = getmodulelist('/SiStrip/MechanicalView/TID/side_2/BadModuleList')
modulesTOB = getmodulelist('/SiStrip/MechanicalView/TOB/BadModuleList')

modulesTECMinusPair = getmodulelistdict('/SiStrip/MechanicalView/TEC/side_1/BadModuleList')
modulesTECPlusPair = getmodulelistdict('/SiStrip/MechanicalView/TEC/side_2/BadModuleList')
modulesTIBPair = getmodulelistdict('/SiStrip/MechanicalView/TIB/BadModuleList')
modulesTIDMinusPair = getmodulelistdict('/SiStrip/MechanicalView/TID/side_1/BadModuleList')
modulesTIDPlusPair = getmodulelistdict('/SiStrip/MechanicalView/TID/side_2/BadModuleList')
modulesTOBPair = getmodulelistdict('/SiStrip/MechanicalView/TOB/BadModuleList')

modulesTECMinusPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TEC/side_1/BadModuleList')
modulesTECPlusPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TEC/side_2/BadModuleList')
modulesTIBPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TIB/BadModuleList')
modulesTIDMinusPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TID/side_1/BadModuleList')
modulesTIDPlusPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TID/side_2/BadModuleList')
modulesTOBPairMap = getmodulereasonsdict('/SiStrip/MechanicalView/TOB/BadModuleList')

#tkmap info
allNBMMap = modulesTECPlusPairMap
allNBMMap.update(modulesTECMinusPairMap)
allNBMMap.update(modulesTIBPairMap)
allNBMMap.update(modulesTIDPlusPairMap)
allNBMMap.update(modulesTIDMinusPairMap)
allNBMMap.update(modulesTOBPairMap)

#recovered modules
modulesTECPlus_recovered = []
modulesTECMinus_recovered = []
modulesTIB_recovered = []
modulesTIDPlus_recovered = []
modulesTIDMinus_recovered = []
modulesTOB_recovered = []

module_filename = "modules_pprun.txt"
if ( RunKey.lower() == 'cosmic_run' ):
    module_filename = "modules_cosmicrun.txt"

ifile_modules = open( str(module_filename) , "r" )

old_tecm_json = json.loads(ifile_modules.readline())
old_tecp_json = json.loads(ifile_modules.readline())
old_tib_json  = json.loads(ifile_modules.readline())
old_tidm_json = json.loads(ifile_modules.readline())
old_tidp_json = json.loads(ifile_modules.readline())
old_tob_json  = json.loads(ifile_modules.readline())

modulesTECMinus_old = []
modulesTECPlus_old = []
modulesTIB_old = []
modulesTIDMinus_old = []
modulesTIDPlus_old = []
modulesTOB_old = []

for mod in old_tecm_json:
    modulesTECMinus_old.append(str(mod))
for mod in old_tecp_json:
    modulesTECPlus_old.append(str(mod))
for mod in old_tib_json:
    modulesTIB_old.append(str(mod))
for mod in old_tidm_json:
    modulesTIDMinus_old.append(str(mod))
for mod in old_tidp_json:
    modulesTIDPlus_old.append(str(mod))
for mod in old_tob_json:
    modulesTOB_old.append(str(mod))

ofile_modules = open( module_filename ,"w")
    
#drop the full lists into files for later comparison
ofile_modules.write(json.dumps(modulesTECMinus)+"\n")
ofile_modules.write(json.dumps(modulesTECPlus)+"\n")
ofile_modules.write(json.dumps(modulesTIB)+"\n")
ofile_modules.write(json.dumps(modulesTIDMinus)+"\n")
ofile_modules.write(json.dumps(modulesTIDPlus)+"\n")
ofile_modules.write(json.dumps(modulesTOB)+"\n")
    
#compare old modules with new list, keep only outliers
for thismodule in modulesTECPlus_old:
    if str(int(thismodule)) in modulesTECPlus:
        modulesTECPlus.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTECPlus_recovered.append(thismodule)

for thismodule in modulesTECMinus_old:
    if str(int(thismodule)) in modulesTECMinus:
        modulesTECMinus.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTECMinus_recovered.append(thismodule)
        
for thismodule in modulesTIB_old:
    if str(int(thismodule)) in modulesTIB:
        modulesTIB.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTIB_recovered.append(thismodule)

for thismodule in modulesTIDPlus_old:
    if str(int(thismodule)) in modulesTIDPlus:
        modulesTIDPlus.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTIDPlus_recovered.append(thismodule)

for thismodule in modulesTIDMinus_old:
    if str(int(thismodule)) in modulesTIDMinus:
        modulesTIDMinus.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTIDMinus_recovered.append(thismodule)

for thismodule in modulesTOB_old:
    if str(int(thismodule)) in modulesTOB:
        modulesTOB.remove( str(int(thismodule)) )
        del allNBMMap[str(int(thismodule))]
    else:
        modulesTOB_recovered.append(thismodule)

n_mTECp = len(modulesTECPlus)
n_mTECm = len(modulesTECMinus)
n_mTIB  = len(modulesTIB)
n_mTIDp = len(modulesTIDPlus)
n_mTIDm = len(modulesTIDMinus)
n_mTOB  = len(modulesTOB)
n_mTotal = n_mTECp + n_mTECm + n_mTIB + n_mTIDp + n_mTIDm + n_mTOB

module_limit = int(newbadmod_limits['newbadmodules'])
isBadModuleProblem = 0
    
if ( ( n_mTECp > module_limit ) or ( n_mTECm > module_limit ) or ( n_mTIB > module_limit ) or ( n_mTIDp > module_limit ) or ( n_mTIDm > module_limit ) or ( n_mTOB > module_limit ) ):
    isProblem = 1
    isBadModuleProblem = 1
    ofile.write("Check: Bad Modules \n")
    ofile.write("\n")

    ofile.write("Number of new bad modules per partition:\n")
    ofile.write("TECPlus: "+str(n_mTECp)+"\n")
    ofile.write("TECMinus: "+str(n_mTECm)+"\n")
    ofile.write("TIB: "+str(n_mTIB)+"\n")
    ofile.write("TIDPlus: "+str(n_mTIDp)+"\n")
    ofile.write("TIDMinus: "+str(n_mTIDm)+"\n")
    ofile.write("TOB: "+str(n_mTOB)+"\n")
#    ofile.write("Link to the current bad-module tracker map: http://vocms01/event_display/TkDocFeedback/newbadmodules_tkmap.png"+"\n\n")
    
    ofile.write("\nBad module details:\n\n")
    if n_mTotal < 1000:
        if ( len(modulesTECPlus) > 0 ):
            ofile.write("TECPlus:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTECPlus:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTECPlusPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")
        if ( len(modulesTECMinus) > 0 ):
            ofile.write("TECMinus:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTECMinus:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTECMinusPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")
        if ( len(modulesTIB) > 0 ):
            ofile.write("TIB:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTIB:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTIBPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")
        if ( len(modulesTIDPlus) > 0 ):
            ofile.write("TIDPlus:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTIDPlus:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTIDPlusPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")
        if ( len(modulesTIDMinus) > 0 ):
            ofile.write("TIDMinus:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTIDMinus:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTIDMinusPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")
        if ( len(modulesTOB) > 0 ):
            ofile.write("TOB:\n")
            ofile.write("--------------------------"+"\n")
            for singlemodule in modulesTOB:
                p = subprocess.Popen("~/scratch0/CMSSW_4_2_3/src/DBINFO/locateModule_mod.sh "+str(singlemodule) , shell=True, stdout=subprocess.PIPE)
                out = p.stdout.read().strip()
                ofile.write(str(singlemodule)+"\n")
                ofile.write("Failure: "+str(modulesTOBPair[singlemodule])+"\n")
                ofile.write(str(out)+"\n")
                ofile.write("--------------------------"+"\n")
            ofile.write("\n")

allRMMap = []
allRMMap += modulesTECPlus_recovered
allRMMap += modulesTECMinus_recovered
allRMMap += modulesTIB_recovered
allRMMap += modulesTIDPlus_recovered
allRMMap += modulesTIDMinus_recovered
allRMMap += modulesTOB_recovered

#dump the new bad modules to a flat file for TrackerMap script
ofileBadM = open("tkmap_bm_input.txt","w")
if isBadModuleProblem == 1:
    for detid, value in allNBMMap.iteritems():
        #mapinfo = str(detid) + " " + str(value)
        mapinfo = str(detid) + " 2"
        ofileBadM.write( mapinfo + "\n" )
    for detid in allRMMap:
        mapinfo = str(detid) + " 1"
        ofileBadM.write( mapinfo + "\n" )
ofileBadM.close()

nr_mTECp = len(modulesTECPlus_recovered)
nr_mTECm = len(modulesTECMinus_recovered)
nr_mTIB  = len(modulesTIB_recovered)
nr_mTIDp = len(modulesTIDPlus_recovered)
nr_mTIDm = len(modulesTIDMinus_recovered)
nr_mTOB  = len(modulesTOB_recovered)
nr_mTotal = nr_mTECp + nr_mTECm + nr_mTIB + nr_mTIDp + nr_mTIDm + nr_mTOB

if ( nr_mTotal > 0 ):
    #not a problem
    ofile.write("\nCheck: Recovered Modules \n")
    ofile.write("\n")
    ofile.write("Number of recovered modules per partition:\n")
    ofile.write("TECPlus: "+str(nr_mTECp)+"\n")
    ofile.write("TECMinus: "+str(nr_mTECm)+"\n")
    ofile.write("TIB: "+str(nr_mTIB)+"\n")
    ofile.write("TIDPlus: "+str(nr_mTIDp)+"\n")
    ofile.write("TIDMinus: "+str(nr_mTIDm)+"\n")
    ofile.write("TOB: "+str(nr_mTOB)+"\n\n")

    if nr_mTotal < 1000:
        ofile.write("Recovered module details:\n\n")
        if ( len(modulesTECPlus_recovered) > 0 ):
            ofile.write("TECPlus:\n")
            for singlemodule in modulesTECPlus_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        if ( len(modulesTECMinus_recovered) > 0 ):
            ofile.write("TECMinus:\n")
            for singlemodule in modulesTECMinus_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        if ( len(modulesTIB_recovered) > 0 ):
            ofile.write("TIB:\n")
            for singlemodule in modulesTIB_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        if ( len(modulesTIDPlus_recovered) > 0 ):
            ofile.write("TIDPlus:\n")
            for singlemodule in modulesTIDPlus_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        if ( len(modulesTIDMinus_recovered) > 0 ):
            ofile.write("TIDMinus:\n")
            for singlemodule in modulesTIDMinus_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        if ( len(modulesTOB_recovered) > 0 ):
            ofile.write("TOB:\n")
            for singlemodule in modulesTOB_recovered:
                ofile.write(str(singlemodule)+"\n")
            ofile.write("\n")
        ofile.write("-------------------------------------------- \n \n")

ofile.close()
ofile_SMS.close()

ofilee = open("sendemail.txt","w")
ofilee.write(str(isProblem)+"\n")
ofilee.close()

if ( isAllFEDsInError == 1 or isNoStripDataProblem == 1 or ( BeamMode.upper() == 'STABLE BEAMS' and IsCollisionsRun == 1 and stripHv == 1 ) ):
    ofilesms = open("sendsms.txt","w")
    ofilesms.write(str(isProblemSMS)+"\n")
    ofilesms.close()
else:
    print "No conditions for SMS"


isFEDProblem = max(wasFEDProblem, isFEDProblem)
isChannelProblem = max(wasChannelProblem, isChannelProblem)
isActiveChannelProblem = max(wasActiveChannelProblem, isActiveChannelProblem)
isAPVeTimingProblem = max(wasAPVeTimingProblem, isAPVeTimingProblem)
isDetFractionProblem = max(wasDetFractionProblem, isDetFractionProblem)
isDetStoNProblem = max(wasDetStoNProblem, isDetStoNProblem)
isTkProblem = max(wasTkProblem, isTkProblem)
isNoStripDataProblem = wasNoStripDataProblem + isNoStripDataProblem

#DEBUG
print "isNoStripDataProblem to write: " , isNoStripDataProblem

ofile_reference = open("reference.txt","w")
ofile_reference.write(str(RunNumber)+"\n")
ofile_reference.write(str(FillNumber)+"\n")
ofile_reference.write(str(thistime)+"\n")
ofile_reference.write(str(isFEDProblem)+"\n")
ofile_reference.write(str(isChannelProblem)+"\n")
ofile_reference.write(str(isActiveChannelProblem)+"\n")
ofile_reference.write(str(isAPVeTimingProblem)+"\n")
ofile_reference.write(str(isDetFractionProblem)+"\n")
ofile_reference.write(str(isDetStoNProblem)+"\n")
ofile_reference.write(str(isTkProblem)+"\n")
ofile_reference.write(str(isNoStripDataProblem)+"\n")
ofile_reference.close()

print "Analysis completed" 

