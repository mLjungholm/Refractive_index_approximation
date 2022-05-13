function flag = yesno(msg)
aceptY = {'yes','Yes','YES','y','Y'};
aceptN = {'No','NO','n','N'};
while true
    prompt = strcat(string(msg)," ",'[Y/N]:');
    awn = input(prompt,'s');
    if any(cellfun(@(x)isequal(x,awn),aceptY))
        flag = true;
        return
    elseif any(cellfun(@(x)isequal(x,awn),aceptN))
        flag = false;
        return
    else
        fprintf('Error: %s is not a valid response \n',string(awn))
    end
end
end