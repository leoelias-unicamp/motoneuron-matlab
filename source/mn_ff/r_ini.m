function synapse = r_ini(indexR,w,time)
%R_INI Fraction of bound receptors of inhibitory synapse.
% SYNAPSE = R_INI(I,W,T) returns a value in order to calculate the synaptic
% current reaching a compartment.

global ti0 s j1 X_ini SYNAPSE_ini t3 pd2 r0_ini aux1
global r_inf_ini tau_R_ini beta_R_ini

if time >= ti0 && time < 0.0
    if s*round((time-w(j1(indexR)))/s) == 0;
        X_ini(indexR) = 1.0;
        j1(indexR) = j1(indexR)+1;
    else
        X_ini(indexR) = 0;
    end

    if X_ini(indexR) == 1
        SYNAPSE_ini(indexR) = 1.0;
        t3(indexR) = time;
        r0_ini(indexR) = aux1(indexR);
    end
    if SYNAPSE_ini(indexR) == 1 && (time-t3(indexR)) > pd2
        SYNAPSE_ini(indexR) = 0.0;
        t3(indexR) = time;
        r0_ini(indexR) = aux1(indexR);
    end

    if SYNAPSE_ini(indexR) > 0
        synapse = r_inf_ini+(r0_ini(indexR)-r_inf_ini)*exp(-(1/tau_R_ini)*(time-t3(indexR)));
    else
        synapse = r0_ini(indexR)*exp(-beta_R_ini*(time-t3(indexR)));
    end
else
    synapse = 0.0;
end

aux1(indexR) = synapse;