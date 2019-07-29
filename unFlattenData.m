function out = unFlattenData(input, x, y)

if size(input,1) ~= length(x) * length(y)
    disp("Error!!")
    out = 3;
else
    oldshape = size(input);
    newshape = [length(x) length(y) oldshape(2:end)];
    out = zeros(newshape);
    for i = 1:length(y)
        out(:,i,:,:,:,:,:,:,:) = input(length(x)*(i-1)+1:length(x)*i,:,:,:,:,:,:,:);
    end
end

end