function L_prime = LEpreprocessor(M,V)

n=length(V);
% the Laplacian Matrix
L=M;
for i=1:n
    d=0;
    for j=1:n
        if M(i,j)~=0 && i~=j
           d = d+1; 
        end
    end
    L(i,:) = (-1)*L(i,:)./d;
    L(i,i)=1;
end

L_prime=[ L zeros(n) zeros(n) 
                zeros(n) L zeros(n)
                zeros(n) zeros(n) L];
            
delta = L*V;


for i=1:n
     ring=[];
     %the neighbors of i
     for   j=1:n
         if M(i,j)~=0
             ring =[ring, j];
         end
     end

    V_i = V(ring,:)';
    
    %the coeff matrix for the system that solves for T
    
    A=zeros(3*length(ring),7);
    for r=1:length(ring)
        A(r,:) = [V_i(1,r) 0 V_i(3,r) (-1)*V_i(2,r) 1 0 0];
        A(r+length(ring),:) = [V_i(2,r) (-1)*V_i(3,r) 0 V_i(1,r) 0 1 0];
        A(r+2*length(ring),:) = [V_i(3,r) V_i(2,r) (-1)*V_i(1,r) 0 0 0 1];
    end
    
    delta_i = delta(i,:)';
    delta_ix = delta_i(1);
    delta_iy = delta_i(2);
    delta_iz = delta_i(3);
    
    Ainv = pinv(A);
    s = Ainv(1,:);
    h = [Ainv(2,:)
            Ainv(3,:)
            Ainv(4,:)];
    t = [Ainv(5,:)
            Ainv(6,:)
            Ainv(7,:)];
    Tdelta = [s*delta_ix-h(3,:)*delta_iy+h(2,:)*delta_iz
                    h(3,:)*delta_ix+s*delta_iy-h(1,:)*delta_iz
                    (-1)*h(2,:)*delta_ix+h(1,:)*delta_iy+s*delta_iz];
    
    %Updating the weights in Lx_prime,Ly_prime,Lw_prime
    L_prime(i,[ring (ring+n) (ring+2*n)]) = (-1)*L_prime(i,[ring (ring+n) (ring+2*n)])+...
                                                                    Tdelta(1,:);
     L_prime(i+n,[ring (ring+n) (ring+2*n)]) = (-1)*L_prime(i+n,[ring (ring+n) (ring+2*n)]) +...
                                                                    Tdelta(2,:);
      L_prime(i+2*n,[ring (ring+n) (ring+2*n)]) =  (-1)*L_prime(i+2*n,[ring (ring+n) (ring+2*n)])+...
                                                                    Tdelta(3,:);
end

end