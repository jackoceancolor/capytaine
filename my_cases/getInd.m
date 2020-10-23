function ind = getInd(dimStruct, str2find)
    ind = 0;
    for j=1:length(dimStruct)
        if string(dimStruct(j).Name) == str2find
            ind = j;
        end
    end
end