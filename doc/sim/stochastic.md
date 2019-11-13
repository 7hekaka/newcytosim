# Stochastic Events in Cytosim
 
Most stochastic events are simulated as such, using pseudo random numbers generated on the fly by [Mersenne Twister](https://en.wikipedia.org/wiki/Mersenne_Twister).

If a stochastic event occurs at a constant rate, its time of occurence in generated to follow an [exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution), using the standard method:

	time = -log(random()) / rate

where `random()` returns a random number uniformly distributed in ]0,1] (zero is excluded). This applies the [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling) method to generate exponential variates. The resulting variable `time` is exponentially distributed with expectancy `1/rate`.


# Bell's Law and Kramers' Reaction Rate Theory

The detachment of a molecular link follows Bell's law. Essentially, the detachment rate `off_rate` varies exponentially with the magnitude of the force experienced by the link:

	off_rate = unbinding_rate * exp( force / unbinding_force )
 
where `force` is the norm of the force vector calculated by `cytosim`, while `unbinding_rate` and `unbinding_force` are constant parameters associated with the bound state. In cytosim, these parameters are specified in the definitions of `Hands`.

> [Models for the specific adhesion of cells to cells](http://dx.doi.org/10.1126/science.347575)  
> Bell, G. I. (1978) - Science, 200(4342), 618–627. 

This law was explained by the theory of Hendrik Kramers:

>  [Brownian motion in a field of force and the diffusion model of chemical reactions ](https://www.sciencedirect.com/science/article/pii/S0031891440900982)  
>  H.A. Kramers - Physica VII, no 4, pp284-304 - 1940

The same law can be expressed differently:

	off_rate = unbinding_rate * exp( force * molecular_scale / ( kB * temperature ) )

Where `kB` is [Boltzman's constant](https://en.wikipedia.org/wiki/Boltzmann_constant).
For more information, see:

> The load dependence of rate constants.  
> Sam Walcott - J Chem Phys - 2008

# Time-varying Rates

Using Kramers rate theory implies that the rate of the event is varying in time.  
The Gillespie approach needs to be modified, and we follow the procedure described in:
 
>  [A Dynamical Monte Carlo Algorithm for Master Equations with Time-Dependent Transition Rates](http://dx.doi.org/10.1007/BF02765541)  
>  A. Prados et al. Journal of Statistical Physics, Vol. 89, Nos. 3/4, 1997  
 
In practice, a normalized time `esp` is first generated,
again using a random number uniformly distributed in [0,1] provided by `random()`.
At each time step, `esp` is reduced as a function of the value of the rate during the interval.
The associated event is performed if `esp` becomes negative, which is when the associated time
crosses the 'present' into the past.
  
	time = 0;
	esp = -log(random());
	
	while ( time < max_time )
	{
		time = time + time_step;
		esp = esp - time_step * rate(time);
		while ( esp < 0 )
		{
		    do_event();
		    esp = esp - log(random())
		}
	}
 
This code must be adapted depending on the circumstances. In this example, the event can be performed multiple times in the same time step, but this would not be done for detachment and other events that can only occur once. 

13.11.2019
