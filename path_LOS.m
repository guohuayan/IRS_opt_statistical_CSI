function [ loss ] = path_LOS( d )
%     d=d/1000;
%     loss=89.5 + 16.9*log10(d);
    loss=38.46 + 20*log10(d);
%     loss=35.3 + 37.6*log10(d);
end

