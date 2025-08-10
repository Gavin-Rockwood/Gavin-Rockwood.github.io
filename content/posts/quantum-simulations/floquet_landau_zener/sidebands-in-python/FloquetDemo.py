import numpy as np
import qutip as qt
import scqubits as scq
import matplotlib.pyplot as plt
import scipy.optimize as opt

import copy
import tqdm
import sys

def state_tracking(state_history, other_sorts = [], reference_states = []):
    """
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
    """


    if len(reference_states) == 0:
        reference_states = state_history[0]

    history = []
    history_other_sorts = []
    previous_states = reference_states
    overlap_history = []
    for state in tqdm.tqdm(range(len(reference_states)), desc = 'Reference States'):
        psi_im1 = reference_states[state]
        history.append([])
        history_other_sorts.append([])
        overlap_history.append([])
        for step in range(len(state_history)):
            overlaps = []
            for j in range(len(state_history[step])):
                overlaps.append(abs(state_history[step][j].dag()*psi_im1))
    
            max_loc = np.argmax(overlaps)
            history[state].append(state_history[step][max_loc])
            to_append = []
            for sort in other_sorts:
                to_append.append(sort[step][max_loc])
            history_other_sorts[state].append(to_append)
            overlap_history[state].append(overlaps)
            psi_im1 = state_history[step][max_loc]

    return {"history": history, "other_sorts": history_other_sorts, "overlap_history": overlap_history}


def Floquet_Sweep(drive_params, ham, drive_op, reference_states = [], sample_times = []):
    """
    Computes Floquet modes and quasienergies for a set of drive parameters, and sorts the results according to reference states.
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
        The sorted Floquet modes and associated quasienergies, as returned by `state_tracking`.
    """

    # Unless sample times are provided, we assume the sampling is at t=0
    if len(sample_times) == 0:
        sample_times = [0]*len(drive_params)
    
    floquet_bases = []
    floquet_modes = []
    floquet_energies = []
    completed_params = {}
    for i in tqdm.tqdm(range(len(drive_params)), desc = "Drive Parameters"):
        # Check if the drive parameters have already been computed
        if str(drive_params[i]) in completed_params.keys():
            floquet_modes.append(floquet_bases[completed_params[str(drive_params[i])]].mode(sample_times[i]))
            floquet_energies.append(floquet_energies[completed_params[str(drive_params[i])]])

        else:
            epsilon = drive_params[i]["epsilon"]
            freq = drive_params[i]["frequency"]
            drive_coef = lambda t: epsilon*np.sin(2*np.pi*freq*t)

            H_T = qt.QobjEvo([ham, [drive_op, drive_coef]])
            floquet_basis = qt.FloquetBasis(2*np.pi*H_T, 1/freq)
            floquet_bases.append(floquet_basis)
            floquet_modes.append(floquet_basis.mode(sample_times[i]))
            floquet_energies.append(floquet_basis.e_quasi)
            completed_params[str(drive_params[i])] = i
    
    sort_res = state_tracking(floquet_modes, other_sorts = [floquet_energies], reference_states = reference_states)

    return sort_res    

def Find_Resonance(ham, drive_op, drive_freqs, epsilon, reference_states, show_plot = True):
    """
    Computes and fits the resonance in a driven quantum system using Floquet theory.
    This function sweeps over a range of Stark shifts, calculates the corresponding quasienergies
    using Floquet analysis, and fits the difference between reference states to a model function.
    Optionally, it plots the quasienergies and fitted difference.
    Parameters
    ----------
    ham : array-like or object
        The system Hamiltonian.
    drive_op : array-like or object
        The operator representing the drive.
    base_freq : float
        The base frequency of the drive (in GHz).
    epsilon : float
        The drive amplitude.
    stark_shifts : array-like
        Array of Stark shift values (in GHz) to sweep over.
    reference_states : list or array-like
        States used as references for extracting quasienergies.
    show_plot : bool, optional
        If True, displays a plot of the quasienergies and fitted difference (default is True).
    Returns
    -------
    list
        [fitted_stark_shift, approximate_drive_time], where:
            fitted_stark_shift : float
                The Stark shift value at the minimum difference (in GHz).
            approximate_drive_time : float
                The approximate drive time (in ns) extracted from the fit.
    Notes
    -----
    - Requires the functions `Floquet_Sweep` and `opt.curve_fit` to be defined/imported.
    - Uses matplotlib for plotting.
    - Assumes quasienergies are returned in `floquet_sweep_res["other_sorts"]`.
    """
    drive_params = []
    for drive_freq in drive_freqs:
        drive_params.append({"frequency": drive_freq, "epsilon": epsilon})
    
    floquet_sweep_res = Floquet_Sweep(drive_params, ham, drive_op, reference_states = reference_states)
    #return floquet_sweep_res

    q_vals1 = np.array([floquet_sweep_res["other_sorts"][0][i][0] for i in range(len(floquet_sweep_res["other_sorts"][0]))])/np.pi
    q_vals2 = np.array([floquet_sweep_res["other_sorts"][1][i][0] for i in range(len(floquet_sweep_res["other_sorts"][1]))])/np.pi
    #q_vals1[q_vals1<0] += 2*(drive_freqs[q_vals1<0])
    #q_vals2[q_vals2<0] += 2*(drive_freqs[q_vals2<0])

    absdifs = abs(q_vals1-q_vals2)

    to_fit = [min([absdifs[i], 2*(drive_freqs[i])-absdifs[i]]) for i in range(len(absdifs))]

    fit_func = lambda x, p0, p1, p2: p2*np.sqrt((x-p0)**2 + p1**2)
    p0 = drive_freqs[np.argmin(to_fit)]
    p1 = np.min(to_fit)
    p2 = abs((absdifs.max() - absdifs.min())/(drive_freqs[np.argmax(to_fit)] - drive_freqs[np.argmin(to_fit)]))
    p = [p0, p1, p2]

    fit_res = opt.curve_fit(fit_func, drive_freqs, to_fit, p0 = p)

    print(f"Fitted Drive Frequency: {fit_res[0][0]} GHz")
    print(f"Approximate Drive Time: {1/(fit_res[0][1]*fit_res[0][2])} ns")
    
    if show_plot:
        fig, ax1 = plt.subplots(figsize=(8, 4), dpi = 200)

        # Plot q_vals1 and q_vals2 on the left y-axis
        ax1.set_xlabel("Drive Frequency (GHz)")
        ax1.set_ylabel(r"Quasienergies/$\pi$ (GHz)")
        ax1.plot(drive_freqs, q_vals1, 'o:', label="State 1", color="forestgreen", lw = 0.5)
        ax1.plot(drive_freqs, q_vals2, 's:', label="State 2", color="dodgerblue", lw = 0.5)
        #ax1.tick_params(axis="y", labelcolor="tab:blue")

        # Create a twin y-axis for difs
        ax2 = ax1.twinx()
        ax2.set_ylabel("Difference (GHz)", color="firebrick")
        ax2.plot(drive_freqs, to_fit, 'd:', label="Difference", color="firebrick", lw = 0.5)

        x = np.linspace(drive_freqs.min(), drive_freqs.max(), 100)
        y = [fit_func(x[i], fit_res[0][0], fit_res[0][1], fit_res[0][2]) for i in range(len(x))]
        ax2.plot(x, y, color="firebrick", lw = 5, alpha = 0.2)
        ax2.tick_params(axis="y", labelcolor="tab:red")
        
        
        ax1.legend(loc="upper right")

        fig.tight_layout()
        plt.title("Stark Shift Fitting")
        plt.show()
    return [fit_res[0][0], 1/(fit_res[0][1]*fit_res[0][2])]


def get_FLZ_flattop(H_op, drive_op, freq, epsilon, envelope_func, ramp_time, psi0, psi1, 
                    num_t_samples=10, epsilons_to_sample=[], n_theta_samples=100, dt=0):
    """
    Python conversion of the Julia get_FLZ_flattop function.
    
    Parameters:
    -----------
    H_op : qutip.Qobj
        The Hamiltonian operator
    drive_op : qutip.Qobj  
        The drive operator
    freq : float
        The drive frequency
    epsilon : float
        The drive amplitude
    envelope_func : function
        The envelope function that takes time as argument
    ramp_time : float
        The total ramp time
    psi0, psi1 : qutip.Qobj
        The two states to track
    num_t_samples : int, optional
        Number of time samples (default: 10)
    epsilons_to_sample : list, optional
        List of epsilon values to sample (if empty, will be generated)
    n_theta_samples : int, optional
        Number of theta samples for optimization (default: 100)
    dt : float, optional
        Time step (if 0, will be set to 1/freq)
        
    Returns:
    --------
    float
        The absolute value of θr/(2π*floq_frequency)
    """
    
    if dt == 0:
        dt = 1/freq
        
    times_to_sample = np.linspace(0, ramp_time, num_t_samples)
    
    if len(epsilons_to_sample) != len(times_to_sample):
        epsilons_to_sample = [[epsilon * envelope_func(t)] for t in times_to_sample]
    
    # Create drive parameters for Floquet sweep
    drive_params = []
    for eps_list in epsilons_to_sample:
        drive_params.append({"epsilon": eps_list[0], "frequency": freq})
    
    # Perform Floquet sweep using existing function
    states_to_track = [psi0, psi1]
    floq_sweep_res = Floquet_Sweep(drive_params, H_op, drive_op, 
                                   reference_states=states_to_track, 
                                   sample_times=times_to_sample)
    
    # Extract quasienergies for the final step
    final_step_idx = len(epsilons_to_sample) - 1
    quasi_energy_0 = floq_sweep_res["other_sorts"][0][final_step_idx][0]
    quasi_energy_1 = floq_sweep_res["other_sorts"][1][final_step_idx][0]
    
    floq_frequency = (quasi_energy_0 - quasi_energy_1) / (2 * np.pi)
    
    # Get final Floquet states
    psi0_floq = floq_sweep_res["history"][0][final_step_idx]
    psi1_floq = floq_sweep_res["history"][1][final_step_idx]
    
    # Create time-dependent Hamiltonian for evolution
    drive_coef = lambda t: epsilon * envelope_func(t) * np.sin(2 * np.pi * freq * t)
    H_drive = qt.QobjEvo([H_op, [drive_op, drive_coef]])
    
    # Solve time evolution
    times_evolution = np.arange(times_to_sample[0], times_to_sample[-1] + dt, dt)
    drive_res_0 = qt.sesolve(2 * np.pi * H_drive, psi0, times_evolution)
    psi0_final = drive_res_0.states[-1]
    
    # Optimization function to minimize
    def to_minimize(theta):
        theta_val = theta[0] if isinstance(theta, (list, np.ndarray)) else theta
        combination = psi0_floq + np.exp(1j * theta_val) * psi1_floq
        overlap = psi0_final.dag() * combination
        return 1 - abs(overlap)**2 / 2
    
    # Grid search for initial guess
    thetas = np.linspace(0, 2*np.pi, n_theta_samples)
    theta_values = [to_minimize(theta) for theta in thetas]
    theta_guess = thetas[np.argmin(theta_values)]
    
    # Optimize
    result = opt.minimize(to_minimize, [theta_guess], method='BFGS')
    theta_opt = result.x[0]
    
    # Calculate θr
    theta_r = np.mod(np.pi - 2 * theta_opt, 2 * np.pi)
    
    print(f"floq_frequency: {floq_frequency}")
    print(f"θ: {theta_opt}, θr: {theta_r}")
    
    return abs(theta_r / (2 * np.pi * floq_frequency))



def Bump_Envelope(t, drive_time, k=2, center=None):
    if center is None:
        center = drive_time / 2
    x = (t - center) / (drive_time / 2)
    if x <= -1 or x >= 1:
        return 0
    elif x == 0:
        return 1
    else:
        return np.exp(k * x**2 / (x**2 - 1))

def Bump_Ramp_Envelope(t, drive_time, ramp_time=1, k=2):
    if t < ramp_time:
        return Bump_Envelope(t, 2 * ramp_time, k=k)
    elif ramp_time <= t <= (drive_time - ramp_time):
        return 1
    elif t > (drive_time - ramp_time):
        return Bump_Envelope(t, 2 * ramp_time, k=k, center=drive_time - ramp_time)
    