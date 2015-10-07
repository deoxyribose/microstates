function microstate_distance = GMD(A,B,nbchan)

if size(A,1) ~= nbchan
    A = A';
    B = B';
end
one = bsxfun(@rdivide,bsxfun(@minus,A,mean(A)),std(A,1)); % standardize
two = bsxfun(@rdivide,bsxfun(@minus,B,mean(B)),std(B,1));
if size(one,2) == 1 % compare one microstate to many
    inde_i_tuborg = bsxfun(@minus,one,two);
elseif size(one,2) == size(two,2) % compare the first microstate in one matrix to the first in the second matrix, etc.
    inde_i_tuborg = one-two;
else
    error('Dimension mismatch')
end
microstate_distance = sqrt(bsxfun(@dot,inde_i_tuborg,inde_i_tuborg)/nbchan);

end