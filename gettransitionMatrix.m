function transitionMatrix = gettransitionMatrix(EEG);
tmp = {};
tmp{1} = EEG.idx;
transitions = cellfun(@(x)([x(1:length(x)-1); x(2:length(x))]), tmp, 'UniformOutput', false);
alltransitions = cell2mat(transitions)';
[uniqueTransitions, ~, i]=unique(alltransitions,'rows','stable');
v=arrayfun(@(x) sum(i==x),1:size(uniqueTransitions,1))';
p = v/sum(v);
transitionMatrix = sparse(uniqueTransitions(:,1), uniqueTransitions(:,2), p,size(uniqueTransitions,1),size(uniqueTransitions,1))
end