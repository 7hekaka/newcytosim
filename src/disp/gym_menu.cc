// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_menu.h"
#include "glut.h"


int gym::createMenu(void (*func)(int))
{
    return glutCreateMenu(func);
}


void gym::addMenuEntry(char const* str, int val)
{
    glutAddMenuEntry(str, val);
}


void gym::addSubMenu(char const* str, int val)
{
    glutAddSubMenu(str, val);
}


void gym::clearMenu(int menu)
{
    glutSetMenu(menu);
    const int x = glutGet(GLUT_MENU_NUM_ITEMS);
    for ( int m = x; m > 0; --m )
        glutRemoveMenuItem(m);
}


void gym::attachMenu(int b)
{
    if ( b==GLUT_LEFT_BUTTON || b==GLUT_RIGHT_BUTTON )
        glutAttachMenu(b);
}

