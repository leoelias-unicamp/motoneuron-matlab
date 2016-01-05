function synapse = r_exc(indexR,w,time)
%R_EXC Fraction of bound receptors of excitatory synapse.
% SYNAPSE = R_EXC(I,W,T) returns a value in order to calculate the synaptic
% current reaching a compartment.

global te0 s j0 X_exc SYNAPSE_exc t2 pd1 r0_exc aux0
global r_inf_exc tau_R_exc beta_R_exc

if j0(indexR) > length(w)
    synapse = 0;
else
    if time >= te0 && time < 0.0
        if s*round((time-w(j0(indexR)))/s) == 0;
            X_exc(indexR) = 1.0;
            j0(indexR) = j0(indexR)+1;
        else
            X_exc(indexR) = 0;
        end
        
        if X_exc(indexR) == 1
            SYNAPSE_exc(indexR) = 1.0;
            t2(indexR) = time;
            r0_exc(indexR) = aux0(indexR);
        end
        if SYNAPSE_exc(indexR) == 1 && (time-t2(indexR)) > pd1
            SYNAPSE_exc(indexR) = 0.0;
            t2(indexR) = time;
            r0_exc(indexR) = aux0(indexR);
        end
        
        if SYNAPSE_exc(indexR) > 0
            synapse = r_inf_exc+(r0_exc(indexR)-r_inf_exc)*exp(-(1/tau_R_exc)*(time-t2(indexR)));
        else
            synapse = r0_exc(indexR)*exp(-beta_R_exc*(time-t2(indexR)));
        end
    else
        synapse = 0.0;
    end
end

aux0(indexR) = synapse;