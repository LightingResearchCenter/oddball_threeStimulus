#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.76.00), Tue 04 Jun 2013 11:10:14 AM EDT
and further modified directly from the created .py file by Petteri Teikari, 2013
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

# PROBLEMS (atm)
# - First run with fresh compile (after changes, fresh .pyc I guess), the sounds are not output

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0

from psychopy import visual, core, data, event, logging, gui, sound, prefs
from psychopy.constants import *  # things like STARTED, FINISHED

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle

import os  # handy system and path functions
import time, sys # for pyo
from datetime import datetime # for timing  if needed

# Use pyo rather than pygame, should have lower latencies
from pyo import * # https://groups.google.com/forum/#!topic/psychopy-users/ZPdusXWPB_A

# EXPERIMENT PARAMETERS
noOfBlockLoops = 1# default is 4 (original article, 6)
noOfIregularTargetLoops = 40 # default is 8


soundVolume = 0.1 # volume, calibrate to 70 db SPL
soundDuration = 1.9925 # empirical value to get ~2000 ms 

triggerDuration = 0.200 # not on during the whole std/odd sound, easier to separate
                                  # successive standard tones when there are gaps between the triggers

introDuration = 2.0 + 5.0 # wav duration + the pause after seconds
introFileWav = 'auditoryStimuli/introSpeak.wav'
outroDuration = 2.0 # seconds, you have to change manually if you use longer WAVs
outroFileWav = 'auditoryStimuli/outroSpeak.wav'

playFromWAVs = 1 # generate "on-fly" (0), or play from WAVs saved on disk (1)
usePyo = 0 # pyo still a bit buggy, and does not provide consistent behavior, with sporadic errors, and some distortion added occasionally to sound

# Function to init reshuffle trials in the irregular target condition
def reShuffleList(trialList):

    # We need to CONSTRAIN HERE, the list so that no 2 oddballs are for example right after each other 
    # i.e. we could "freeze" certain elements of the list, specified in the input .csv -file
    # e.g. http://stackoverflow.com/questions/12238005/python-shuffling-list-but-keeping-some-elements-frozen
    
    # a bit non-elegant way to extract to the flags for the freeze from the csv-file
    '''
    for index, item in enumerate(irregularTargets.trialList):        
        print "index = ", index
        freezeBooleanFlag = item.freeze        
        print "freeze = ", freezeBooleanFlag    
    '''
     
    # memorize position of fixed elements   
    fixed = [(pos, item) for (pos,item) in enumerate(trialList) if item.freeze]
    # print " "
    # print "FIXED: \n", fixed
    # print " "
    
    # shuffle list
    random.shuffle(trialList)
    # print "Shuffle, 1st Pass: \n", trialList
    # print " "
    
    # swap fixed elements back to their original position
    for pos, item in fixed:        
        index = trialList.index(item)
        trialList[pos], trialList[index] = trialList[index], trialList[pos]
    
    # print "Shuffled List \n", trialList # print shuffled irregularTargets
    # print " "
    return trialList
    
        # NOTE! (comment columns added to make lines unique)
        # http://stackoverflow.com/questions/12238005/python-shuffling-list-but-keeping-some-elements-frozen
        # @Pawel: this solution might break if there are duplicates in the list. Also (less important) items.index() can scan the whole list for a single value. It makes the algorithm quadratic in time. – J.F. Sebastian Sep 2 '12 at 18:20
        # @J.F. Sebastian You are right, thanks for pointing this out! However, for the problem at hand (shuffling answers to quiz questions) neither should pose a problem. – tobias_k Sep 2 '12 at 18:32
        
        # NOTE2!: This shuffling with constraints actually only addresses the problem with consecutive oddballs, but does not
        # allow as long as 8-12 consecutive standard tones without oddballs

# Function to init PYO if used
def initPyo():

    # CREATE THE pyo SERVER
    pyo = Server(sr=48000, nchnls=1, buffersize=512, duplex=0) # , audio='jack', jackname='pyo')
        # Latency is `buffer size / sampling rate` in seconds. Defaults to 256.

    # print debug info
    # print pa_get_output_devices() # prints the detected output devices
        # 'Music Streamer II: USB Audio (hw:1,0)'
        # 'Music Streamer II: USB Audio (hw:2,0)', 'sysdefault', 'front', 'rear', 'surround40', 'iec958', 'spdif', 'pulse', 'dmix', 'default'], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
    # print " "
    # print pa_get_devices_infos()
        # 'name': 'Music Streamer II: USB Audio (hw:1,0)'}, 3: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100
        # 'name': 'Music Streamer II: USB Audio (hw:2,0)'}, 6: {'host api index': 0, 'latency': 0.04265306144952774, 'default sr': 44100, 'name': 'sysdefault'}, 7: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'front'}, 8: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'rear'}, 9: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'surround40'}, 10: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'iec958'}, 11: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'spdif'}, 12: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'pulse'}, 13: {'host api index': 0, 'latency': 0.04266666620969772, 'default sr': 48000, 'name': 'dmix'}, 14: {'host api index': 0, 'latency': 0.01160997711122036, 'default sr': 44100, 'name': 'default'}})
    # print " "
    # print "Recognized devices: \n", pa_list_devices()
        # 11: IN, name: spdif, host api index: 0, default sr: 44100 Hz, latency: 0.011610 s
        # 11: OUT, name: spdif, host api index: 0, default sr: 44100 Hz, latency: 0.011610 s

    # Define parameters 
    print " "
    pyo.setAmp(soundVolume) # http://www.iact.umontreal.ca/pyo/manual/Server.html
    # pyo.setOutputDevice(14) # Set the audio output device number. See `pa_list_devices()`
    # pyo.setInOutDevice(5)
    # pa_count_devices()

        # for setting the default output device for ALSA, see:
        # https://wiki.archlinux.org/index.php/ALSA#Set_the_default_sound_card
            # sudo gedit /etc/modprobe.d/alsa-base.conf        
                # specify the order of sound cards:
                # options snd slots=snd_usb_audio,snd_cmipci,snd_hda_intel
                # options snd_usb_audio index=0 # USB-Audio - Music Streamer II www.hirestech.com 2010 REV 1.7 Music Streamer II at usb-0000:00:1d.7-5.1, full
                # options snd_cmipci index=1 # CMI8738 - C-Media CMI8738 (model 37) at 0xdc00, irq 16
                # options snd_hda_intel index=2 # HDA-Intel at 0xfe9dc000 irq 45

    # pyo.setMidiInputDevice(5) #: Set the MIDI input device number. See `pm_list_devices()`.
    # pyo.setMidiOutputDevice(5) #: Set the MIDI output device number. See `pm_list_devices()`.
        # Might be more elegant to check for the name spdif and Music Streamer II somehow
        # rather than hard-coding this number here

    print " "
    print "Initializing pyo (for playing audio) ..."
    print "     (If you bunch of warning see messages below or above, they are not critical.)"
    print "     see e.g. http://blog.yjl.im/2012/11/pyaudio-portaudio-and-alsa-messages.html for details"
    print "End of PyAudio initialization."
    pyo.boot()

        # for an example and some discussion see for example:
        # https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw

    # from: https://groups.google.com/forum/#!topic/psychopy-users/qNKAHIVPG2I
    # Following 3 lines apparently do not affect to anything, override by File -> Preferences
        # psychopy.prefs.general['audioLib'] = ['pygame'] 
        # psychopy.prefs.general['audioDriver'] = ['portaudio']     
        # print(psychopy.prefs) #tell me about all current settings 

    print " "
    print "Server created: ", serverCreated()
    print "Default output (PortAudio): ", pa_get_default_output() # Returns the index number of Portaudio's default output device.
    print "Default output (PortMidi): ", pm_get_default_output() # Returns the index number of Portmidi's default output device.
    print " "
    print pa_list_devices()
    print " "

    # To see the manual of pyo, visit: 
    # http://www.iact.umontreal.ca/pyo/manual/
        
    # To get the pyo running at Ubuntu: compile it from the source (http://code.google.com/p/pyo/wiki/Installation)
        # sudo python setup.py install --use-jack
        # you may need to install some development packages:
        # sudo apt-get install libportmidi-dev libsndfile1-dev portaudio19-dev 
            # To install pyo, you will need the following dependencies:
            # Python 2.6.x or 2.7.x
            # portaudio - portaudio19-dev 
            # portmidi - libportmidi-dev
            # libsndfile - libsndfile1-dev
            # liblo - Install Manually (http://liblo.sourceforge.net/), and read INSTALL
            # Ubuntu: python-tk (doesn't come installed by default on Ubuntu)
            # Ubuntu: python-dev package

    # Possible errors: "Portaudio error: Invalid number of channels"
        # check: https://groups.google.com/forum/#!topic/psychopy-users/LoFzw7F3yx0
        # edit the sound.py from "/usr/local/lib/python2.7/dist-packages/PsychoPy-1.77.00-py2.7.egg/psychopy" (or where PsychoPy is installed)
        # sudo gedit /usr/local/lib/python2.7/dist-packages/PsychoPy-1.77.00-py2.7.egg/psychopy/sound.py

    # To test the JACK/pyo functioning, run "QjackCtl" from "Sound & Video"
        # you may need to add the line "@audio -       rtprio          99" to
            # sudo gedit /etc/security/limits.conf
        # and add yourself 
            # sudo usermod -a -G audio <userName>

    # For some resources to optimize the real-time performance of your Linux System, see:
        # http://wiki.linuxaudio.org/wiki/system_configuration#hpet
        # http://www.alsa-project.org/main/index.php/Low_latency_howto
        # http://wiki.linuxaudio.org/wiki/real_time_info
        
        # Tuning up for USB Audio Card
            # http://alsa.opensrc.org/Usb-audio

    # For enabling/disabling services on your Linux, check out:
    # http://askubuntu.com/questions/19320/recommended-way-to-enable-disable-services

        # sudo apt-get install sysv-rc-conf
        # sudo sysv-rc-conf
        
        # or GUI-based
        # sudo apt-get install jobs-admin
        
        # or BUM (http://askubuntu.com/questions/131684/how-to-boot-with-bluetooth-turned-off)
        # sudo apt-get install bum && sudo bum
        
        # See also:
        # http://www.hecticgeek.com/2012/06/few-things-to-speed-up-ubuntu/

    # You might get a bunch of ALSA lib.xxx error messages in Ubuntu, while the sound works ok,
    # to solve them, see for example the following thread:
    # http://stackoverflow.com/questions/7088672/pyaudio-working-but-spits-out-error-messages-each-time
    # sudo gedit /usr/share/alsa/alsa.conf
    
        # cat /proc/asound/cards
        # 2 [II             ]: USB-Audio - Music Streamer II
        #                  www.hirestech.com 2010 REV 1.7 Music Streamer II at usb-0000:00:1d.7-5.1, full 

        # aplay  -l
        # **** List of PLAYBACK Hardware Devices ****
        # card 2: II [Music Streamer II], device 0: USB Audio [USB Audio]
        #   Subdevices: 1/1
        #   Subdevice #0: subdevice #0

    # If you get the following error with maybe Segmentation Fault
    # ALSA lib pcm_dmix.c:1018:(snd_pcm_dmix_open) unable to open slave
    # check the following: http://chakra-project.org/bbs/viewtopic.php?pid=50909
    # sudo gedit ~/.asoundrc (II from cat /proc/asound/cards, remember to boot after modifying
        # pcm.!default {
        # type hw
        # card II
        # }

        # ctl.!default {
        # type hw
        # card II
        # }    

    # Init the sound stimulus, and update only in the loop
    # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw
    pyo.start()
    return pyo


if usePyo == 1:
    pyo = initPyo()        
else:
    print "pygame used instead of 'pyo'" # for pygame

# Import LabJack U6 USB Daq
import u6 

#initalize LabJack U6
d  = u6.U6(debug = False) # add try/catch here
print "ConfigU6: ", d.configU6()

# Configure the LabJack U6
print "ConfigIO: ", d.configIO()
print " "

# DEFINE THE REGISTERS, if needed
# not needed with our .getFeedback() -call
    # FI01 for IRREGULAR Boolean
    # FI02 and FI03 for TRIGGERS (Standard and Deviant)
    # http://labjack.com/support/modbus/ud-modbus
    # http://labjack.com/support/labjackpython

    # Init button values (LabJack now not reading buttons, sent directly to the 
    # Biosemi ActiveTwo Trigger input.
    # d.getFeedback(u6.BitStateWrite( 0, 0 ), u6.BitStateWrite( 1, 0 ) )
    # d.getFeedback(u6.BitDirWrite( 0, 0 ), u6.BitDirWrite( 1, 0 ) )

# Store info about the experiment session
expName = 'Oddball_initial'  # from the Builder filename that created this script
expInfo = {u'session': u'001', u'participant': u''}
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

    
# Setup files for saving
if not os.path.isdir('data'):
    os.makedirs('data')  # if this fails (e.g. permissions) we will get error
filename = 'data' + os.path.sep + '%s_%s' %(expInfo['participant'], expInfo['date'])
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)

# Setup the Window
win = visual.Window(size=(280, 150), fullscr=False, screen=0, allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-1.000,-1.000,-1.000], colorSpace='rgb')

# Initialize components for Routine "Introduction"
IntroductionClock = core.Clock()

if usePyo == 1:
    Intro_Speak = SfPlayer(introFileWav ,speed=1,loop=False)  # http://www.iact.umontreal.ca/pyo/manual/SfPlayer.html
else:
    Intro_Speak = sound.Sound(u'A', secs=1.0) # only works with pygame
    Intro_Speak.setVolume(soundVolume) # only works with pygame

intro_textField = visual.TextStim(win=win, ori=0, name='intro_textField',
    text=u'Intro',    font=u'Arial',
    pos=[0, 0], height=0.3, wrapWidth=None,
    color=u'white', colorSpace=u'rgb', opacity=1,
    depth=-1.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()

if usePyo == 1:
    # Unfinished feature, finish at some point and see if there is an effect for
    # latency and jitter
    if playFromWAVs == 1:
        soundStimulus = SfPlayer('auditoryStimuli/standardTone.wav',speed=1,loop=False)  # http://www.iact.umontreal.ca/pyo/manual/SfPlayer.html
    else:
        soundStimulus = snd = Sine(freq=1000, mul=soundVolume)
else:
    soundStimulus = sound.Sound('A', secs=1)
    soundStimulus.setVolume(soundVolume)
    
# Initialize components for Routine "Outro"
OutroClock = core.Clock()

if usePyo == 1:
    Outro_Speak= SfPlayer(outroFileWav,speed=1,loop=False)  # http://www.iact.umontreal.ca/pyo/manual/SfPlayer.html
else:
    Outro_Speak = sound.Sound(u'A', secs=1.0)  # only works with pygame
    Outro_Speak.setVolume(soundVolume)  # only works with pygame    
    
outro_textField = visual.TextStim(win=win, ori=0, name='outro_textField',
    text=u'End of experiment',    font=u'Arial',
    pos=[0, 0], height=0.3, wrapWidth=None,
    color=u'white', colorSpace=u'rgb', opacity=1,
    depth=-1.0)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# Define the Cycle loop for THREE-STIMULUS ODDBALL PARADIGM

#------Prepare to start Routine "Introduction"-------
t = 0
IntroductionClock.reset()  # clock 
frameN = -1
routineTimer.add(1.000000)
# update component parameters for each repeat
# keep track of which components have finished
IntroductionComponents = []
IntroductionComponents.append(Intro_Speak)
IntroductionComponents.append(intro_textField)
for thisComponent in IntroductionComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "Introduction"-------
continueRoutine = True
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = IntroductionClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # start/stop Intro_Speak
    if usePyo == 1:
        Intro_Speak.out()
        if frameN >= 0 and Intro_Speak.isPlaying == False:
            # keep track of start time/frame for later
            Intro_Speak.tStart = t  # underestimates by a little under one frame
            Intro_Speak.frameNStart = frameN  # exact frame index
            Intro_Speak.play()  # start the sound (it finishes automatically)
        elif Intro_Speak.isPlaying == True and t >= (Intro_Speak.tStart + introDuration):                
            Intro_Speak.stop()  # stop the sound (if longer than duration) 
    else:
        Intro_Speak.setSound(introFileWav)
        if frameN >= 0 and Intro_Speak.status == NOT_STARTED:
            # keep track of start time/frame for later
            Intro_Speak.tStart = t  # underestimates by a little under one frame
            Intro_Speak.frameNStart = frameN  # exact frame index
            Intro_Speak.play()  # start the sound (it finishes automatically)
        elif Intro_Speak.status == STARTED and t >= (Intro_Speak.tStart + introDuration):
            Intro_Speak.stop()  # stop the sound (if longer than duration)        
    
    
    # *intro_textField* updates
    if t >= 0.0 and intro_textField.status == NOT_STARTED:
        # keep track of start time/frame for later
        intro_textField.tStart = t  # underestimates by a little under one frame
        intro_textField.frameNStart = frameN  # exact frame index
        intro_textField.setAutoDraw(True)
    elif intro_textField.status == STARTED and t >= (0.0 + introDuration):
        intro_textField.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in IntroductionComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the [Esc] key)
    if event.getKeys(["escape"]):
        core.quit()

    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

    d.getFeedback(u6.BitStateWrite( 0, 1 ) ) # FI00 (0), i.e. START RECORD Trigger
    core.wait(introDuration, hogCPUperiod=0)  # need hogCPUperiod < wait time 
                                        # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw

#-------Ending Routine "Introduction"-------
for thisComponent in IntroductionComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)


# START THE BLOCK
CycleLoop = data.TrialHandler(nReps=noOfBlockLoops, method=u'sequential', 
    extraInfo=expInfo, originPath=None,
    trialList=[None],
    seed=None, name='CycleLoop')

startTime = datetime.now() # Get start time for SOA/ISI saving
        # probably more elegant to update at some point to use the handy timers created above
timePrev = startTime # init the timing value

thisExp.addLoop(CycleLoop)  # add the loop to the experiment
thisCycleLoop = CycleLoop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb=thisCycleLoop.rgb)

if thisCycleLoop != None:
    for paramName in thisCycleLoop.keys():
        exec(paramName + '= thisCycleLoop.' + paramName)

print " "
print "NEW CYCLE"      
    
for thisCycleLoop in CycleLoop:
    currentLoop = CycleLoop
    # print "LOOP COUNT (Cycle No): ", currentLoop
    # abbreviate parameter names if possible (e.g. rgb = thisCycleLoop.rgb)
    if thisCycleLoop != None:
        for paramName in thisCycleLoop.keys():
            exec(paramName + '= thisCycleLoop.' + paramName)
    
    print " "
    print "NEW LOOP (Irregular)"       
    
    # LOOP for IRREGULAR TARGETS
    irregularTargets = data.TrialHandler(nReps=noOfIregularTargetLoops, method=u'sequential', # NOTE! Randomized after this "manually" with constraints
        extraInfo=expInfo, originPath=None,
        trialList=data.importConditions('Oddballs_threeStimulus_WavFiles.csv'),
        seed=None, name='irregularTargets')

    irregularTargets.data.addDataType('ISI')#this will help store ISI with the stimuli, taken from demo/TrialHandler.py
    d.getFeedback(u6.BitStateWrite( 1, 1 )) # Write irregular boolean is HIGH for FO01       

        # NOTE now how there are only 4 repetitions (nReps) of the irregularTargets,
        # whereas you might assume 8 to be the correct number, but now as we define
        # the deviant-deviant interval differently so that the deviant can be preceded by
        # 2-6 or 8-12 standard tones and during 8 consecutive tones there might be two deviants
        # for example the each block now consists of 16 tones in contrast to 8 tones per block
        # of the regular condition (see below), see the Fig. 1 of Jongsma et al. (2013) for graphical
        # representation        
        
    # print "Original List \n", irregularTargets.trialList # print irregularTargets
    # print "Length of list: ", len(irregularTargets.trialList)
    # print type(irregularTargets.trialList)
    
    thisExp.addLoop(irregularTargets)  # add the loop to the experiment
    thisIrregularTarget = irregularTargets.trialList[0]  # so we can initialise stimuli with some values
    
    #print thisExp
        # <psychopy.data.ExperimentHandler object at 0xa5d664c>
    # print thisIrregularTarget
        #{'Duration': 0.20000000000000001, 'SOA': 0.59399999999999997, 'oddBallYes': 0, 'Subfolder': u'auditoryStimuli', 'Filename': u'standardTone.wav'}

    # abbreviate parameter names if possible (e.g. rgb=thisIrregularTarget.rgb)
    if thisIrregularTarget != None:
        for paramName in thisIrregularTarget.keys():
            exec(paramName + '= thisIrregularTarget.' + paramName)        
                 
    tonesPerBlock = int(irregularTargets.nTotal / irregularTargets.nReps) # i.e. 8 
    
    for thisIrregularTarget in irregularTargets: # will run nReps * nTrials times (i.e. 6 * 8 = 48 times)
        
        # now we want shuffle the list on each repetition (on first index, =0)        
        if irregularTargets.thisIndex == 0: # make a function out of reshuffling at some point
            # print irregularTargets
            irregularTargets.trialList = reShuffleList(irregularTargets.trialList)
        
        # abbreviate parameter names if possible (e.g. rgb = thisIrregularTarget.rgb)        
        if thisIrregularTarget != None:
            for paramName in thisIrregularTarget.keys():
                exec(paramName + '= thisIrregularTarget.' + paramName)
        
        #------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock 
        frameN = -1
        timeNow= datetime.now()
        # update component parameters for each repeat

        Sounds = os.path.join(Subfolder, Filename) # so that relative path works both in Linux/Mac and Win
        soundStimulus.setVolume(soundVolume)
        # https://groups.google.com/forum/?fromgroups#!topic/psychopy-users/twtPKicPGeQ        
        # print Sounds
        # print thisIrregularTarget
        
        if oddBallYes== 1 and distracterYes == 0:
            
            # what to do when oddball is presented, i.e. send a digital output as a trigger to the EEG system            
            soundStimulus.setSound(Sounds)
            if usePyo == 1:
                soundStimulus.out()
                
            # see, e.g. http://labjack.com//support/labjackpython
            d.getFeedback(u6.BitStateWrite( 3, 1 ), u6.BitStateWrite( 2, 0 ), u6.BitStateWrite( 1, 0) )
            #buttons = d.getFeedback(u6.BitStateRead( 0 ), u6.BitStateRead( 1 ) )
            
            # directions  = d.getFeedback(u6.BitDirRead( 0 ), u6.BitDirRead( 1 ) )
            # ainValue = d.getAIN(0)  # Read from AIN0 in one function    
            # print "      button 1 (FI00): ", buttons[0], "        direction: ", directions[0]
            # print ",     button 2 (FI01): ", buttons[1], "        direction: ", directions[1]
            # print "      AIN00: ", "%1.2f" % ainValue, " V"            
            
            ISI_ms = 1000*(timeNow-timePrev).total_seconds()
            # print "      Irregular oddball, ISI:", ISI_ms, " ms"
            
            print "      Irregular oddball, ISI:", ISI_ms, " ms", ", index: ", irregularTargets.thisIndex+1, "/", tonesPerBlock, ", N: ", irregularTargets.thisN, "/", irregularTargets.nTotal
            irregularTargets.data.add('ISI', ISI_ms)
            
            if usePyo == 1:
                # time.sleep() # SOA: column in the CSV-file
                core.wait(Duration+SOA, hogCPUperiod=0)  # need hogCPUperiod < wait time 
                                        # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw
                                        # optimally would have "+ .getDurationOfWav or something"
            
        elif oddBallYes==0 and distracterYes == 0:

            # what to do when standard tone is presented, i.e. other digital output line trigger to EEG
            soundStimulus.setSound(Sounds)
            if usePyo == 1:                
                soundStimulus.out()
                
            # see, e.g. http://labjack.com//support/labjackpython
            d.getFeedback(u6.BitStateWrite( 3, 0 ), u6.BitStateWrite( 2, 1 ), u6.BitStateWrite( 1, 0) )
            # buttons = d.getFeedback(u6.BitStateRead( 0 ), u6.BitStateRead( 1 ) )
            # directions  = d.getFeedback(u6.BitDirRead( 0 ), u6.BitDirRead( 1 ) )
            # ainValue = d.getAIN(0)  # Read from AIN0 in one function    
            # print "      button 1 (FI00): ", buttons[0], "        direction: ", directions[0]
            # print ",     button 2 (FI01): ", buttons[1], "        direction: ", directions[1] 
            # print "      AIN00: ", "%1.2f" % ainValue, " V"                      

            ISI_ms = 1000*(timeNow-timePrev).total_seconds()
            print "   Irregular standard, ISI:", ISI_ms, " ms"
            irregularTargets.data.add('ISI', ISI_ms)
            
            if usePyo == 1:
                # time.sleep() # SOA: column in the CSV-file
                core.wait(Duration+SOA, hogCPUperiod=0)  # need hogCPUperiod < wait time 
                                        # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw
                                        # optimally would have "+ .getDurationOfWav or something"
        
        elif distracterYes == 1:
            
            # what to do when standard tone is presented, i.e. other digital output line trigger to EEG
            soundStimulus.setSound(Sounds)
            if usePyo == 1:                
                soundStimulus.out()
                
            # see, e.g. http://labjack.com//support/labjackpython
            d.getFeedback(u6.BitStateWrite( 3, 0 ), u6.BitStateWrite( 2, 0 ), u6.BitStateWrite( 1, 1) )
            # buttons = d.getFeedback(u6.BitStateRead( 0 ), u6.BitStateRead( 1 ) )
            # directions  = d.getFeedback(u6.BitDirRead( 0 ), u6.BitDirRead( 1 ) )
            # ainValue = d.getAIN(0)  # Read from AIN0 in one function    
            # print "      button 1 (FI00): ", buttons[0], "        direction: ", directions[0]
            # print ",     button 2 (FI01): ", buttons[1], "        direction: ", directions[1] 
            # print "      AIN00: ", "%1.2f" % ainValue, " V"                      

            ISI_ms = 1000*(timeNow-timePrev).total_seconds()
            print "           Irregular distracter, ISI:", ISI_ms, " ms", ", index: ", irregularTargets.thisIndex+1, "/", tonesPerBlock, ", N: ", irregularTargets.thisN, "/", irregularTargets.nTotal
            irregularTargets.data.add('ISI', ISI_ms)
            
            if usePyo == 1:
                # time.sleep() # SOA: column in the CSV-file
                core.wait(Duration+SOA, hogCPUperiod=0)  # need hogCPUperiod < wait time 
                                        # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw
                                        # optimally would have "+ .getDurationOfWav or something"
            
            
        else:
            print "   Definition of sound incorrect for irregular targets!"
        
        # keep track of which components have finished
        trialComponents = []
        # trialComponents.append(LearningOddballParadigm) # for pygame
        trialComponents.append(soundStimulus)
        
        # print "Volume: ", soundStimulus.getVolume()
        
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        timePrev = timeNow

        #-------Start Routine "trial"-------
        continueRoutine = True
        triggerOff = 0
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            if usePyo == 1:
                if frameN >= 0 and soundStimulus.isPlaying == False:
                    # keep track of start time/frame for later
                    soundStimulus.tStart = t  # underestimates by a little under one frame
                    soundStimulus.frameNStart = frameN  # exact frame index
                    soundStimulus.play()  # start the sound (it finishes automatically)
                # elif soundStimulus.isPlaying == True and  t >= triggerDuration:
                #    d.getFeedback(u6.BitStateWrite( 3, 0 ), u6.BitStateWrite( 2, 0 ) ) # std and odd LOW
                elif soundStimulus.isPlaying == True and triggerOff == 0 and t >= (soundStimulus.tStart + soundDuration):                
                    soundStimulus.stop()  # stop the sound (if longer than duration)            
                    triggerOff = 1; # reset flag, only goes once inside this
            else:
                if frameN >= 0 and soundStimulus.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    soundStimulus.tStart = t  # underestimates by a little under one frame
                    soundStimulus.frameNStart = frameN  # exact frame index
                    soundStimulus.play()  # start the sound (it finishes automatically)
                elif soundStimulus.status == STARTED and triggerOff == 0 and  t >= triggerDuration:
                    d.getFeedback(u6.BitStateWrite( 3, 0 ), u6.BitStateWrite( 2, 0 ), u6.BitStateWrite( 1, 0 ) ) # std and odd LOW                    
                    triggerOff = 1; # reset flag, only goes once inside this
                elif soundStimulus.status == STARTED and t >= (soundStimulus.tStart + soundDuration):
                    soundStimulus.stop()  # stop the sound (if longer than duration)
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested that we end
                routineTimer.reset()  # this is the new t0 for non-slip Routines
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # check for quit (the [Esc] key)
            if event.getKeys(["escape"]):
                core.quit()
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()            
        
        #-------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        
        thisExp.nextEntry()
        
    # completed all the repeats of 'irregularTargets'
    # irregularTargets.printAsText(stimOut=[], dataOut=('ISI_mean', 'ISI_std', 'ISI_raw'))
    
    # get names of stimulus parameters
    if irregularTargets.trialList in ([], [None], None):  params = []
    else:  params = irregularTargets.trialList[0].keys()
    # save data for this loop
    irregularTargets.saveAsExcel(filename + '.xlsx', sheetName='irregularTargets',
        stimOut=params,
        dataOut=['n','all_mean','all_std', 'all_raw'])
        

 # completed x repeats of 'CycleLoop'

# get names of stimulus parameters
if CycleLoop.trialList in ([], [None], None):  params = []
else:  params = CycleLoop.trialList[0].keys()
# save data for this loop
CycleLoop.saveAsExcel(filename + '.xlsx', sheetName='CycleLoop',
    stimOut=params,
    dataOut=['n','all_mean','all_std', 'all_raw'])

# END OF THE EXPERIMENT
# PUT ALL THE TRIGGERS to LOW
d.getFeedback(u6.BitStateWrite( 3, 0 ) ) # FI03 (3), set LOW (0), i.e. trigger for ODDBALL 
d.getFeedback(u6.BitStateWrite( 2, 0 ) ) # FI02 (2), set LOW (0), i.e. trigger for STANDARD TONE
d.getFeedback(u6.BitStateWrite( 1, 0 ) ) # FI01 (1), set LOW (0), i.e. trigger for DISTRACTER
d.getFeedback(u6.BitStateWrite( 0, 0 ) ) # FI00 (0), set LOW (0), i.e. START RECORD Trigger

# Print debug info on console (or stats of the experiment)
print " "
print "Experiment took:", (datetime.now() - startTime), " in total"
durationSeconds = (datetime.now() - startTime).total_seconds()
print "Experiment took:", durationSeconds, " seconds in total"
theorDuration = noOfBlockLoops * ((noOfIregularTargetLoops*8))*2.0
print "Theoretical duration: ", theorDuration, " seconds" # number of cycles * (irreg + reg) * SOA
print 100*(durationSeconds / theorDuration), "% of theoretical duration"
print " "

time.sleep(0.5)

#------Prepare to start Routine "Outro"-------
t = 0
OutroClock.reset()  # clock 
frameN = -1
routineTimer.reset()  # if we abort early the non-slip timer needs reset
routineTimer.add(1.000000)
# update component parameters for each repeat
# keep track of which components have finished
OutroComponents = []
OutroComponents.append(Outro_Speak)
OutroComponents.append(outro_textField)

for thisComponent in OutroComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED        

#-------Start Routine "Outro"-------
continueRoutine = True
# print routineTimer.getTime()
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = OutroClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # start/stop Outro_Speak
    if usePyo == 1:
        Outro_Speak.out()
        if frameN >= 0 and Outro_Speak.isPlaying == False:
            # keep track of start time/frame for later
            Outro_Speak.tStart = t  # underestimates by a little under one frame
            Outro_Speak.frameNStart = frameN  # exact frame index
            Outro_Speak.play()  # start the sound (it finishes automatically)
        elif Outro_Speak.isPlaying == True and t >= (Outro_Speak.tStart + OutroDuration):                
            Outro_Speak.stop()  # stop the sound (if longer than duration) 
    else:
        Outro_Speak.setSound(outroFileWav)
        if frameN >= 0 and Outro_Speak.status == NOT_STARTED:
            # keep track of start time/frame for later
            Outro_Speak.tStart = t  # underestimates by a little under one frame
            Outro_Speak.frameNStart = frameN  # exact frame index
            Outro_Speak.play()  # start the sound (it finishes automatically)
        elif Outro_Speak.status == STARTED and t >= (Outro_Speak.tStart + OutroDuration):
            Outro_Speak.stop()  # stop the sound (if longer than duration)      
    
    
    # *outro_textField* updates
    if t >= 0.0 and outro_textField.status == NOT_STARTED:        
        # keep track of start time/frame for later
        outro_textField.tStart = t  # underestimates by a little under one frame
        outro_textField.frameNStart = frameN  # exact frame index
        outro_textField.setAutoDraw(True)
    elif outro_textField.status == STARTED and t >= (0.0 + outroDuration):
        outro_textField.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in OutroComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the [Esc] key)
    if event.getKeys(["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

    core.wait(outroDuration, hogCPUperiod=0)  # need hogCPUperiod < wait time 
                                        # from: https://groups.google.com/forum/#!topic/psychopy-dev/GcObydzJgyw

#-------Ending Routine "Outro"-------
for thisComponent in OutroComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

d.close() # Close the Labjack U6 device

if usePyo == 1:
    pyo.stop() # Stop the audio callback loop
    # pyo.shutdown() #  Shut down and clear the server. This method will erase all objects
            # from the callback loop. This method need to be called before changing 
            # server's parameters like `samplingrate`, `buffersize`, `nchnls`, ...

win.close()
core.quit()

