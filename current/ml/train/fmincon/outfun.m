function stop = outfun(x,optimValues,state)

persistent history
if isempty(history)
  history.x = [];
  history.fval = [];
end

stop = false;

switch state
  case 'init'
  case 'iter'
    % Concatenate current point and objective function
    % value with history. x must be a row vector.
    history.fval = [history.fval; optimValues.fval];
    history.x = [history.x; x];
    
%     fprintf(['\nCost:' num2str(optimValues.fval) '  XP: ' num2str(x) '\n'])
  case 'done'        
  otherwise
end
end