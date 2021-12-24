% This functions is depricated.
%
% Because this function does not print but displays the design, it is now
% called display_design(). Otherwise nothing changed.
%
% For convinience the function calls display_design(cfg) for you.

function table = print_design(cfg)

warning('print_design:deprecated', 'print_design() is deprecated because the name is not correct - it does not print the design, but displays it.\nIt still works, but please use display_design() instead.')

if nargout == 0
    display_design(cfg);
else
    table = display_design(cfg);
end