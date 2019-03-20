mac:
	g++ -Wall -o cytomaster cytomaster.cpp -D__MACOSX_CORE__ RtMidi.cpp -framework CoreMidi -framework CoreFoundation -framework CoreAudio

alsa:
	g++ -Wall -D__LINUX_ALSA__ -o cytomaster cytomaster.cpp RtMidi.cpp -lasound -lpthread 


install: mac
	cp cytomaster ../live/cytomaster
