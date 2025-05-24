// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef PLAYER_H
#define PLAYER_H

///this turns on display code in some objects
#define DISPLAY 1

#include "gym_color.h"
#include "gym_color_list.h"

#include "simul.h"
#include "parser.h"
#include "display.h"
#include "sim_thread.h"
#include "display_prop.h"
#include "player_prop.h"
#include "property_list.h"

class FiberDisp;
class View;

/// Player is a base for Cytosim's 'play' executable
class Player
{
public:

    /// Simulation object
    Simul        simul;
    
    /// the display parameters
    DisplayProp  disp;
    
    /// the parameters for play
    PlayerProp   prop;
    
    /// container for fiber's display properties
    PropertyList mFiberDisp;
    
    /// container for object's display properties
    PropertyList mPointDisp;

    /// SimThread to control the live simulation
    SimThread    worker;
    
    /// the current Display object
    Display  *   mDisplay;
    
    /// a flag for live simulation
    bool goLive;

    //---------------------------------COMMANDS---------------------------------
    
    /// constructor
    Player();

    ///
    ~Player();
    
    /// initialize display
    void initialize();
  
    /// release memory
    void destroy();

    /// ask to redraw the scene
    void refresh();

    /// return all Fiber FiberDisp
    PropertyList allFiberDisp() const;
   
    /// return all Fiber FiberDisp
    PropertyList allVisibleFiberDisp() const;
 
    /// return all Hand PointDisp
    PropertyList allHandDisp() const;
    
    /// return all Hand PointDisp for which 'visible==true'
    PropertyList allVisibleHandDisp() const;

    /// return all Sphere/Solid/Bead PointDisp
    PropertyList allSphereDisp() const;
 
    /// return all Sphere/Solid/Bead PointDisp for which 'visible==true'
    PropertyList allVisibleSphereDisp() const;

    /// return all Space PointDisp
    PropertyList allSpaceDisp() const;
    
    /// return a FiberDisp
    FiberDisp * firstFiberDisp() const;

    //---------------------------------COMMANDS---------------------------------
    
    /// reset view, without changing the current frame
    void rewind();
   
    /// start animation
    bool startPlayback();
    
    /// accelerate animation if possible
    void accelerate();
    
    /// change replay speed
    void setTimelapse(unsigned);
    
    /// start reverse animation
    bool startBackward();
    
    /// start live simulation
    void extendLive();

    /// stop animation
    void stop();

    /// start or stop animation
    void startstop();
    
    /// reset the sim-state and timer
    void restart(int);

    /// load previous frame
    void previousFrame();
    
    /// go to the next frame, returns 1 if EOF is reached
    void nextFrame();
    
    /// write global display parameters
    void writePlayParameters(std::ostream& out, bool prune) const;

    /// write Object display parameters
    void writeDisplayParameters(std::ostream& out, bool prune) const;
    
    //-----------------------------DISPLAY--------------------------------------
  
    /// initialize display with given style
    void setStyle(unsigned);

    /// build message that appears on top
    std::string buildReport(std::string) const;
    
    /// build message that appears on bottom
    std::string buildLabel() const;
    
    /// build central message
    std::string buildMemo(int) const;

    
    /// set View::focus and rotation to track objects in the simulation
    void autoFocus(View&, Simul const&, unsigned) const;
    
    /// adjust the model view and load frame if asked
    void setPixelSize(View&);

    /// adjust the model view and load frame if asked
    void prepareDisplay(View&);
    
    /// read parameters contained in string
    void readDisplayString(View&, std::string const&);
    
    /// draw all simulation components
    void drawCytosim();
    
    /// draw system calling drawCytosim
    void drawSystem(View&);
    
    /// export current viewport to image file 'filename'
    int saveView(View const&, const char* filename, const char* format, int downsample) const;

    /// export current viewport to an image file
    int saveView(View&, size_t indx, int downsample) const;
    
    /// save high-resolution image of the current scene to 'filename'
    int saveScene(int mag, const char* filename, const char* format, int downsample);
    
    /// save high-resolution image of the current scene
    int saveScene(int mag, const char* root, unsigned indx, int downsample=1);

};


#endif
