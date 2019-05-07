// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 @file
 @brief Global Compile Switches
*/

#ifndef SIM_H
#define SIM_H

/**
 Enable code to be able to read old trajectory files
 Option normally ON
 */
#define BACKWARD_COMPATIBILITY


/**
 If the keyword below is defined, the viscous drag of the fibers
 will be different in the transverse and parallel directions, such that
 it will be 2x easier to move a fiber along it longitudinal direction.
 
 This is unpublished development, and you should set to zero
 */
#define NEW_ANISOTROPIC_FIBER_DRAG 0


/**
 Enables Myosin, Kinesin and Dynein
 Option normally OFF
 */
#define NEW_HANDS 0


/** 
 Enables advanced Space
 */
#define NEW_SPACES 1


#endif
