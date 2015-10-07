function vars = get_component_energies(A,S)
    J = size(A,1);
    T = size(S,2);
    Avar=diag(A'*A)/J;
    Svar=diag(S*S')/T;
    vars=Avar.*Svar;
    vars=vars/sum(vars);
end