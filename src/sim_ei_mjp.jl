function sim_ei_mjp(E_init::Int64, I_init::Int64, alpha::Float64, nu::Float64, gamma::Float64, stop_time::Float64)
    # Setting up initial states and frames
    t = 0.0
    num_exposed = E_init
    num_infectious = I_init
    rt = alpha/nu 
    state_frame = [t num_exposed num_infectious alpha rt]
    # simulate till end point reached
    while (num_exposed > 0 || num_infectious > 0) && t < stop_time
        # Time to the next event
        next_event = rand(Exponential(1 / (alpha * num_infectious + gamma * num_exposed + nu * num_infectious)))

        # Update time
        t += next_event

        # Choose which event happens proportional to the rates
        which_event = sample(["infection", "infectious", "recovery"], Weights([alpha * num_infectious, gamma * num_exposed, nu * num_infectious]))

        # If the event is infection
        if which_event == "infection"
            num_exposed += 1
        # If the event is becoming infectious
        elseif which_event == "infectious"  
            # Update states
            num_exposed -= 1
            num_infectious += 1
        # If the event is recovery
        else
            # Update states
            num_infectious -= 1
        end
        # update state_frame
        rt = alpha/nu
        state_frame = vcat(state_frame, [t num_exposed num_infectious alpha rt])
    end
    return (DataFrame(state_frame, [:time, :E, :I, :alpha, :rt]))
end


