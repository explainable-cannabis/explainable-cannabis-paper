%code to make a correlation graph object in matlab
%
%inputs need to be
%1) a cell array of ROIs
%2) fMRI time-series (parcels x samples)
%3) Thresh = number of connections in the lower triangular matrix
function [CORRgraph] = ELR_make_correlation_graph(ROIs,fMRIdat,thresh)

nodenames = ROIs;

func_mat = tril(corr(fMRIdat'),-1);
func_thresh = func_mat>=min(maxk(func_mat(:),thresh));

s={};
t={};

for i = 1:numel(ROIs)
    for j = i:numel(ROIs)
        if func_thresh(j,i)==1
            s{end+1}=nodenames{i};
            t{end+1}=nodenames{j};
        end
    end
end
CORRgraph = graph(s,t,[],nodenames);

end