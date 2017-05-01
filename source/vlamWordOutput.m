% get handles of all open figures
h = get(0,'Children');

% activate each figure
 for i=1:length(h)
     figure(h(i));
     grid on
 end
 