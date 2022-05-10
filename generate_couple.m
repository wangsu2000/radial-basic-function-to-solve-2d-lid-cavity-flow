function C = generate_couple(A,B)%生成张量积对应的序偶
C = [reshape([kron(A',ones(size(B)))],1,length(A)*length(B));reshape([kron(ones(size(A))',B)],length(A)*length(B),1)';];