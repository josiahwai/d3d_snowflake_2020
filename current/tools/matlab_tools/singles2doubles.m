workspacevars = whos;
for i = 1:length(workspacevars)  
  if strcmp(workspacevars(i).class,'single')   
    name = workspacevars(i).name;    
    assignin('base', name, double(evalin('base', name)));    
  end  
end