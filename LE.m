%L_prime (3n,3n)     the matrix returned by preprocessor
%V (n,3)                    original vertex coord
%static_anchors       index of static points
%handle_anchors     index of handle points
%handle_new_pos    new position of handle points
%w                             weight of position coord influence against Laplacian coord

function V_prime= LE(L_prime,V,static_anchors,handle_anchors,handle_new_pos,w)

anchors = [static_anchors handle_anchors];
n=length(V);
V(handle_anchors,:)=handle_new_pos;
U= reshape(V,3*n,1);

anchor_vec = sparse(1,3*n);
anchor_vec([anchors,anchors+n,anchors+2*n]) = 1;
anchor_mat = diag(anchor_vec);

B_prime = anchor_mat'*anchor_mat;
rhs = B_prime*U;


A_prime=L_prime+B_prime;

V_final = A_prime\rhs;
V_prime = [V_final(1:n) V_final((n+1):(2*n)) V_final((2*n+1):(3*n)) ];
end