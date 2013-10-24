% L_prime is the Matrix returned by Preprocessor
function V_final_xyz = LE(L_prime,V,static_anchors,handle_anchors,handle_anchors_pos)
anchors = [static_anchors handle_anchors];
V_pi = V;
n=length(V);
V_pi(handle_anchors,:) = handle_anchors_pos;
A_prime=L_prime;
rhs = zeros(3*n,1);
for j=1:length(anchors)
    A_prime=[A_prime
        (1:(3*n))==anchors(j)
        (1:(3*n))==(anchors(j)+n)
        (1:(3*n))==(anchors(j)+2*n) ];
    rhs=[rhs
        V_pi(anchors(j),1)
        V_pi(anchors(j),2)
        V_pi(anchors(j),3)];
end

V_final = A_prime\rhs;
rhs_re = A_prime*V_final;
V_final_xyz = [V_final(1:n) V_final((n+1):(2*n)) V_final((2*n+1):(3*n)) ];
end