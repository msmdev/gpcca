function [ optval, opt ] = getopt( opt, field, default )
% get or set structure field values
    %
    % [ optval, opt ] = getopt( opt, field, default )
    %
    % Input:
    %   opt         structure
    %   field       structure field name string
    %   default     field value to set, if the field doesnt exist yet
    % Output:
    %   optval      variable with value equal the field value
    %   opt         (updated) structure 
    % from NLSCON v2.3.3 by Nowak and Weimann
    
    if isfield(opt,field)
        optval = getfield(opt,field) ;
    else
        optval = default ;
        opt = setfield(opt,field,default) ;
    end
end

