    This was run with:
    Python version: 3.11.10 (main, Nov 26 2024, 19:12:45) [Clang 16.0.0 (clang-1600.0.26.4)]
    Qutip version: 5.1.1
    scqubits version: 4.3


## Some Initial Tools
Before looking at Floquet-Landau-Zener dynamics, we need a function designed to sweep floquet modes for a list of parameters. 
### Step 1: State Tracking
In order to understand our system, it will be important to understand how the states evolve as we change parameters in our system. In particular, we will use it to look at the evolution of the Floquet quasienergies and search for avoided crossings! We do this with the function `state_tracking`

    
        Tracks the evolution of quantum states for some list of states at different steps (general sweeping parameters).
        Args:
            state_history (list of list): A list where each element is a list of quantum states at a given time step.
            other_sorts (list of list, optional): Additional lists of data associated with each state, to be tracked alongside the main state history.
            reference_states (list, optional): List of reference quantum states to track. If empty, uses the first set of states in `state_history`.
        Returns:
            dict: A dictionary containing:
                - "history": List of lists of tracked quantum states for each reference state.
                - "other_sorts": List of lists of associated data for each tracked state.
                - "overlap_history": List of lists of overlap values between reference states and all states at each time step.
        


### Step 2: Floquet Sweep
To sweep Floquet modes (using `FloquetBasis` from `qutip` for the for the floquet modes), we will use the function `floquet_sweep`. This function will take a list of parameters and return the Floquet modes for each parameter.

    
        Computes Floquet modes and quasi-energies for a set of drive parameters, and sorts the results according to reference states.
        Parameters
        ----------
        drive_params : list of dict
            List of dictionaries specifying drive parameters for each sweep point. Each dictionary should contain at least
            'epsilon' (drive amplitude) and 'frequency' (drive frequency).
        ham : qutip.Qobj
            The static part of the system Hamiltonian.
        drive_op : qutip.Qobj
            The operator representing the driven part of the Hamiltonian.
        reference_states : list, optional
            List of reference states used for sorting the Floquet modes. Default is an empty list.
        sample_times : list, optional
            List of times at which to sample the Floquet modes for each drive parameter. If not provided, defaults to t=0 for all.
        Returns
        -------
        sort_res : tuple
            The sorted Floquet modes and associated quasi-energies, as returned by `state_tracking`.
        


## Finding Resonances
To find resonances, we will use the function `Find_Resonance`. This function will sweep over a range of Stark shifts, calculate the corresponding quasienergies using Floquet analysis, and fit the difference between reference states to a model function. 

The avoided crossing will have the form
$$
\Delta E(\nu) = a\sqrt{(\nu  - b)^2 + c}
$$
where $\Delta E(\nu)$ is the difference in quasienergies, $\nu$ is the drive frequency, and $a$, $b$, and $c$ are fitting parameters.
