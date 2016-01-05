function synapse = r_ini(i,w,t)
%R_INI Fraction of bound receptors of inhibitory synapse.
% SYNAPSE = R_INI(I,W,T) returns a value in order to calculate the synaptic
% current reaching a compartment.

global ti0 di s aux2 j1 X_ini SYNAPSE_ini t3 pd2 r0_ini aux3
global r_inf_ini tau_R_ini beta_R_ini

if t >= ti0 && t < 0.0
    if s*round((t-(ti0+aux2(i)+w(j1(i))))/s) == 0;
        X_ini(i) = 1.0;
        aux2(i) = aux2(i)+w(j1(i));
        j1(i) = j1(i)+1;
    else
        X_ini(i) = 0;
    end

    if X_ini(i) == 1
        SYNAPSE_ini(i) = 1.0;
        t3(i) = t;
        r0_ini(i) = aux3(i);
    end
    if SYNAPSE_ini(i) == 1 && (t-t3(i)) > pd2
        SYNAPSE_ini(i) = 0.0;
        t3(i) = t;
        r0_ini(i) = aux3(i);
    end

    if SYNAPSE_ini(i) > 0
        synapse = r_inf_ini+(r0_ini(i)-r_inf_ini)*exp(-(1/tau_R_ini)*(t-t3(i)));
    else
        synapse = r0_ini(i)*exp(-beta_R_ini*(t-t3(i)));
    end
else
    synapse = 0.0;
end

aux3(i) = synapse;