# blochSim
Simple Bloch Equation Simulation for MRI. I made this code after following along with the tutorial given at http://mrsrl.stanford.edu/~brian/bloch/. If you have any questions about how it works, that is a great resource.

Also, I don't give any guarantees that this will work in every situation. Please feel free to report bugs or submit additional methods.
## Standard Available Methods

### addSpin(M0,T1,T2,df,reps)
> Adds spins to initialized blochSim object
### tip(alpha, phase)
> Tip all spins by flip angle alpha (radians) (with excitation phase phase (radians))
### freeprecess(time_step)
> Simulates progression of time (time_step, in seconds) accounting for T1,T2 and df.
### sig = signal()
> Returns the Transverse Signal (using the most recent excitation phase as the reciever phase)
### spoil(low_ext,up_ext)
> Gradient Spoiling. (reps from addSpin must be >0 for this to have any effect)

> default low_ext = -2*pi, default up_ext=2*pi
### forcespoil()
> Cheating. Forces the transverse magnetization to 0;
## Standard Helper Methods

### sig = transverseSignal()
> Returns the Transverse Signal (using the most recent excitation phase as the reciever phase)
### sig = longitudinalSignal()
> Returns the longitudinal Signal
### record(bool_optional)
> Adds current object state as data to object memory 
### save(memory_key)
> Appends current transverse signal, longitudinal signal, and time to object array accesible by obj.mem(memory_key)
### plot_memory(memkeys_to_plot (vector: see save function), plot_longitudinal (default false), plot_phase (defualt false, plots magnitude))
> Plots saved signal over time.
### Mplot()
> Plots M as a 3D vector
### reset()
> Resets the object (ie, clear memory, reset internal variables, freeprecess for a long long time)


