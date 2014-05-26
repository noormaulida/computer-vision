function [result,c]=classification(network, new_inputs, new_targets, targets)
    % Test Network
    outputs = sim(network,new_inputs);
    [c,cm] = confusion(new_targets,outputs);
    
    if ((cm(1,1)==1 || cm(1,2)==1))
        result='Benign';
    elseif (cm(2,1)==1 || cm(2,2)==1)
        result='Malignant';
    else
        result='Unclassified';
    end
end