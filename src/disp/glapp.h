// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
//F. Nedelec, Dec 2005

#ifndef GLAPP_H
#define GLAPP_H

#include "real.h"
#include "vector3.h"
#include "quaternion.h"
#include "view.h"
#include <vector>

///glApp extends GLUT to manipulate a 2D or 3D display window
namespace glApp
{
    /// different View associated with the display windows
    extern std::vector<View> views;
    
    /// state of mouse button: up/down
    extern int mouseS;
    
    /// flag indicating that display should be redone
    extern int drawRefresh;
    
    /// initialize first view
    void initialize();
    
    /// change dimensionnality (this affects mostly the mouse controls availability)
    void setDimensionality(int d);
    
    /// create new window
    int newWindow(const char name[]);

    /// create new window with display functions `func` and `mag`
    int newWindow(const char name[], int (*func)(View&), void (*mag)(View&));

    /// create new window with display function `func`
    int newWindow(const char name[], int (*func)(View&));

    /// create new window with display function `func`
    inline int newWindow(int (*func)(View&)) { return newWindow("Cytosim", func); }

    /// destroy window
    void deleteWindow(int win);
    
    /// enter or exit full-screen mode
    void toggleFullScreen();

    /// maximize window size within the current screen
    void toggleWindowSize();
    
    /// callback function for window resize event
    void resizeWindow(int, int);
    
    /// set the range normally visible for zoom = 1
    void setScale(float);

    /// return view associated with current window
    View& currentView();
    
    /// reset current view
    void resetView();
    
    /// draw System
    void displayOtherWindows(int (*drawFunc)(View&));

    /// save higher resolution image with magnification 'mag'
    int saveImage(const char* name, unsigned mag, unsigned downsample);
    
    //--------------------------------- MENUS -----------------------------------
    
    /// build menu, and attach it if argument is a valid button (-1: do not attach)
    int buildMenu();
    
    /// callback function for menu events
    void processMenuEvent(int item);
    
    /// attach default menu to button
    void attachMenu(int menu = 0);

    //-----------------------------------KEYS------------------------------------

    /// returns a string describing mouse and keyboard commands
    void help(std::ostream&);
    
    /// callback function for normal keys
    void processNormalKey(unsigned char, int modifiers);
    
    /// callback function for normal keys
    void processNormalKey(unsigned char, int mouse_x, int mouse_y);

    /// set callback for keyboard events
    void normalKeyFunc(void (*func)(unsigned char, int, int));
    
    /// callback function for normal keys
    void processSpecialKey(int key, int modifiers);

    /// callback function for arrow/function keys
    void processSpecialKey(int key, int mouse_x, int mouse_y);
    
    /// set callback for special key pressed events
    void specialKeyFunc(void (*func)(int, int, int));

    //-----------------------------------MOUSE-----------------------------------
        
    /// callback function for mouse button down/up
    void processMouseClick(int button, int state, int x, int y);
    
    /// callback function for mouse motion, when a button is pressed
    void processMouseDrag(int x, int y);
    
    /// callback function for mouse motion, when no button is pressed
    void processPassiveMouseMotion(int x, int y);
    
    /// set callback for shift-click, with unprojected down-position
    //func(mouseX, mouseY, mouseDown, specialKeys);
    void actionFunc(void (*func)(int, int, const Vector3&, int));
    
    /// set callback for shift-click, with unprojected down- and current- mouse positions
    //func(mouseX, mouseY, mouseDown, mousePosition, specialKeys);
    void actionFunc(void (*func)(int, int, Vector3&, const Vector3&, int));
    
    //---------------------------------MESSAGES---------------------------------
    
    /// set message displayed on current window
    void setMessage(std::string const&);
    
    /// set message displayed on current window
    void setSubtitle(std::string const&);

    /// display given text on screen for 3 sec
    void flashText(std::string const&);
    
    /// display text for 3 sec on screen, used to signify something to user
    template < typename Arg1, typename... Args >
    void flashText(const char* fmt, Arg1&& arg1, Args&&... args)
    {
        char str[1024] = { 0 };
        snprintf(str, sizeof(str), fmt, arg1, args...);
        flashText(str);
    }

    //-------------------------------DISPLAY------------------------------------
    
    /// called after display of scene
    void displayInfo(int W, int H);

    /// display function for main window
    void displayMain();
    
    /// display function for secondary windows
    void displayPlain();
    
    /// declare that current display needs a refresh
    void postRedisplay();
}


#endif
