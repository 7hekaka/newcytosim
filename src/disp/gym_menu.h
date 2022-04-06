// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University


namespace gym
{
    /// shortcut
    extern bool has_menus;
    
    /// shortcut
    int createMenu(void (*func)(int));
    
    /// shortcut
    void addMenuEntry(char const* str, int val);
    
    /// shortcut
    void addSubMenu(char const* str, int val);
    
    /// shortcut
    void clearMenu(int menu);
    
    /// shortcut
    void attachMenu(int b);

}


