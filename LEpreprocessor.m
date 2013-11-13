%M  adjacency matrix
%V (n,3) vertex matrix
function L_prime = LEpreprocessor(M,V)

n=length(V);

% the Laplacian Matrix
D = 1./double(sum(M));
L=speye(n)-diag(D)*double(M);

[I,J, K]  = find(L);
I = [I I+n I+2*n];
J = [J J+n J+2*n];
L_prime = sparse(I,J,repmat(K,3,1));

delta = full(L)*V;


for i=1:n
    
     %the neighbors of i
    ring=[find(M(i,:)),i];
    n_i = length(ring);
    V_i = V(ring,:)';
    
    %the coeff matrix for the system that solves for T
    
    A=zeros(3*length(ring),7);
    for r=1:length(ring)
        A(r,:) = [V_i(1,r) 0 V_i(3,r) (-1)*V_i(2,r) 1 0 0];
        A(r+n_i,:) = [V_i(2,r) (-1)*V_i(3,r) 0 V_i(1,r) 0 1 0];
        A(r+2*n_i,:) = [V_i(3,r) V_i(2,r) (-1)*V_i(1,r) 0 0 0 1];
    end
    
    delta_i = delta(i,:)';
    delta_ix = delta_i(1);
    delta_iy = delta_i(2);
    delta_iz = delta_i(3);
    
    Ainv = pinv(A);
    s = Ainv(1,:);

    h1 = Ainv(2,:);
    h2 = Ainv(3,:);
    h3 = Ainv(4,:);
    t = [Ainv(5,:)
           Ainv(6,:)
           Ainv(7,:)];
    Tdelta = [s*delta_ix-h3*delta_iy+h2*delta_iz
                    h3*delta_ix+s*delta_iy-h1*delta_iz
                    (-1)*h2*delta_ix+h1*delta_iy+s*delta_iz];


    
    %Updating the weights in Lx_prime,Ly_prime,Lz_prime
    rN2N = [ring ring+n ring+2*n];
    L_prime(i,rN2N) = (-1)*L_prime(i,rN2N) + Tdelta(1,:);
     L_prime(i+n,rN2N) = (-1)*L_prime(i+n,rN2N) + Tdelta(2,:);
      L_prime(i+2*n,rN2N) =  (-1)*L_prime(i+2*n,rN2N) + Tdelta(3,:);
end
L_prime=L_prime'*L_prime;

end