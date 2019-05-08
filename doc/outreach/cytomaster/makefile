mac:
	g++ -Wall -o cytomaster cytomaster.cpp -D__MACOSX_CORE__ RtMidi.cpp rtmidi_c.cpp -framework CoreMidi -framework CoreFoundation -framework CoreAudio

alsa:
	g++ -Wall -D__LINUX_ALSA__ -o cytomaster cytomaster.cpp RtMidi.cpp -lasound -lpthread 


install: mac
	cp cytomaster live/cytomaster

ready:
	cp live/0-live.command live/1-live.command
	cp live/0-live.command live/2-live.command
	cp live/0-live.command live/3-live.command
	cp live/0-live.command live/4-live.command
	cp live/0-live.command live/5-live.command
	cp live/0-live.command live/6-live.command
	cp live/0-live.command live/7-live.command
	cp live/0-live.command live/8-live.command
