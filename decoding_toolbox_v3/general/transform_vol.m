% function newvol = transform_vol(vol)
%
% Transforms a 3D volume to a number of axial 2D slices (i.e. along
% dimensions 1 and 2)
%
% Martin Hebart (2009)
%

function newvol = transform_vol(vol)

if ndims(vol)~=3
    error('Input must be 3D!')
end

sz = size(vol);

factor = sz(2)/sz(1);
nRows = round(sqrt(sz(3))*factor);
nColumns = ceil(sz(3)/nRows);

newvol = zeros(sz(1)*nRows,sz(2)*nColumns);

counterRow = 0;
counterColumn = 1;

slicetemplate = newvol;
slicetemplate(1:numel(slicetemplate)) = 1:numel(slicetemplate);
sliceindices = slicetemplate(1:sz(1),1:sz(2))-1;

for i = 1:sz(3)
    if counterRow == nRows
        counterRow = 0;
        counterColumn = counterColumn + 1;
    end
    counterRow = counterRow + 1;
    currindices = sliceindices+slicetemplate((counterRow-1)*sz(1)+1,(counterColumn-1)*sz(2)+1);
    newvol(currindices) = vol(:,:,i);
end