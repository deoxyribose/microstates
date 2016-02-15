function run_length_encoding_length = description_length(Z)
    if min(size(Z)) > 1
        Z = sum(bsxfun(@times,Z,2.^[0:size(Z,2)-1]),2); % convert binary to decimal
    else
        Z = Z';
    end
    J=find(diff([0; Z]));
    encoding = [Z(J), diff([J; size(Z,1)+1])];
    run_length_encoding_length = size(encoding,1);
end
