function synapse = r_exc(i,w,t)
%R_EXC Fraction of bound receptors of excitatory synapse.
% SYNAPSE = R_EXC(I,W,T) returns a value in order to calculate the synaptic
% current reaching a compartment.

global te0 de s aux0 j0 X_exc SYNAPSE_exc t2 pd1 r0_exc aux1
global r_inf_exc tau_R_exc beta_R_exc

if t >= te0 && t < 5000.0
    if s*round((t-(te0+aux0(i)+w(j0(i))))/s) == 0;
        X_exc(i) = 1.0;
        aux0(i) = aux0(i)+w(j0(i));
        j0(i) = j0(i)+1;
    else
        X_exc(i) = 0;
    end

    if X_exc(i) == 1
        SYNAPSE_exc(i) = 1.0;
        t2(i) = t;
        r0_exc(i) = aux1(i);
    end
    if SYNAPSE_exc(i) == 1 && (t-t2(i)) > pd1
        SYNAPSE_exc(i) = 0.0;
        t2(i) = t;
        r0_exc(i) = aux1(i);
    end

    if SYNAPSE_exc(i) > 0
        synapse = r_inf_exc+(r0_exc(i)-r_inf_exc)*exp(-(1/tau_R_exc)*(t-t2(i)));
    else
        synapse = r0_exc(i)*exp(-beta_R_exc*(t-t2(i)));
    end
else
    synapse = 0.0;
end

aux1(i) = synapse;